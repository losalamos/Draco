//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:07 2000
 * \brief  CDI class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "CDI.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <cmath>

namespace rtt_cdi
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS AND DESTRUCTORS
//---------------------------------------------------------------------------//

CDI::CDI(const std_string &id)
    : matID(id),
      grayOpacities(constants::num_Models, 
		    SF_GrayOpacity(constants::num_Reactions)),
      multigroupOpacities(constants::num_Models,
			  SF_MultigroupOpacity(constants::num_Reactions))
{
    Ensure (grayOpacities.size() == constants::num_Models);
    Ensure (multigroupOpacities.size() == constants::num_Models);
}

//---------------------------------------------------------------------------//
    
CDI::~CDI() 
{
    // empty
}


//---------------------------------------------------------------------------//
// STATIC DATA
//---------------------------------------------------------------------------//

std::vector<double> CDI::frequencyGroupBoundaries = std::vector<double>();

//---------------------------------------------------------------------------//
// STATIC FUNCTIONS
//---------------------------------------------------------------------------//

std::vector<double> CDI::getFrequencyGroupBoundaries()
{
    return frequencyGroupBoundaries;
}

//---------------------------------------------------------------------------//

int CDI::getNumberFrequencyGroups()
{
    int ng = 0;

    if (!frequencyGroupBoundaries.empty())
	ng = frequencyGroupBoundaries.size() - 1;

    return ng;
}    

//---------------------------------------------------------------------------//

// nested unnamed namespace that holds data and services used by the
// Planckian integration routines.  The data in this namespace is accessible
// by the methods in this file only (internal linkage).

namespace
{

// constants used in the Taylor series expansion of the Planckian.
const double coeff_3  =    1.0 / 3.0;
const double coeff_4  =   -1.0 / 8.0;
const double coeff_5  =    1.0 / 60.0;
const double coeff_7  =   -1.0 / 5040.0;
const double coeff_9  =    1.0 / 272160.0;
const double coeff_11 =   -1.0 / 13305600.0;
const double coeff_13 =    1.0 / 622702080.0;
const double coeff_15 =  -6.91 / 196151155200.0;
const double coeff_17 =    1.0 / 1270312243200.0;
const double coeff_19 = -3.617 / 202741834014720.0;
const double coeff_21 = 43.867 / 107290978560589824.0;

const double pi       =    2.0 * std::asin(1.0);
const double coeff    =   15.0 / (pi*pi*pi*pi);

// return the 21-term Taylor series expansion for the normalized Planck
// integral given x
inline double taylor_series_planck(double x)
{
    Require (x >= 0.0);

    // calculate the 21-term Taylor series expansion for x
    double xsqrd  = x * x;
    double xpower = x * x * x;
    double taylor = 0.0;
    {
	taylor += coeff_3  * xpower;
	xpower *= x;
	taylor += coeff_4  * xpower;
	xpower *= x;
	taylor += coeff_5  * xpower;
	xpower *= xsqrd;
	taylor += coeff_7  * xpower;
	xpower *= xsqrd;
	taylor += coeff_9  * xpower;
	xpower *= xsqrd;
	taylor += coeff_11 * xpower;
	xpower *= xsqrd;
	taylor += coeff_13 * xpower;
	xpower *= xsqrd;
	taylor += coeff_15 * xpower;
	xpower *= xsqrd;
	taylor += coeff_17 * xpower;
	xpower *= xsqrd;
	taylor += coeff_19 * xpower;
	xpower *= xsqrd;
	taylor += coeff_21 * xpower;
	
	taylor *= coeff;
    }

    Ensure (taylor >= 0.0);
    return taylor;
}

// return the 10-term Polylogarithmic expansion (minus one) for the Planck
// integral given x
inline double polylog_series_minus_one_planck(double x)
{
    Require (x >= 0.0);

    double xsqrd = x * x;

    // calculate the 10-term Polylogirithmic expansion (minus one) for x
    double poly = 0.0;
    
    // calculate the polylogarithmic terms for x bound
    double li1 = 0.0;
    double li2 = 0.0; 
    double li3 = 0.0;
    double li4 = 0.0;

    for (int i = 1; i <= 10; i++)
    {
	li1 += std::exp(-i * x) / i;
	li2 += std::exp(-i * x) / (i*i);
	li3 += std::exp(-i * x) / (i*i*i);
	li4 += std::exp(-i * x) / (i*i*i*i);
    }

    // calculate the lower polylogarithmic integral
    poly = -1.0 * coeff * (xsqrd * x * li1   +
			   3.0 * xsqrd * li2 + 
			   6.0 * x * li3     +
			   6.0 * li4 );
    
    Ensure (poly <= 0.0);
    return poly;
}

}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Planckian spectrum over a frequency range.
 *
 * This function integrates the normalized Plankian that is defined:
 *
 * b(x) = 15/pi^4 * x^3 / (e^x - 1)
 *
 * where 
 *
 * x = (h nu) / (kT)
 * 
 * and 
 * 
 * B(nu, T)dnu = acT^4/4pi * b(x)dx
 * 
 * where B(nu, T) is the Plankian and is defined
 *
 * B(nu, T) = 2hnu^3/c^2 (e^hnu/kt - 1)^-1
 *
 * The normalized Plankian, integrated from 0 to infinity, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Planck function may integrate to something less than one.
 *
 * This function performs the following integration:
 *
 *      Int_(x_low)^(x_high) b(x) dx
 *
 * where x_low is calculated from the input low frequency bound and x_high is
 * calculated from the input high frequency bound.  This integration uses the
 * method of B. Clark (JCP (70)/2, 1987).  We use a 10-term Polylogarithmic
 * expansion for the normalized Planckian, except in the low-x limit, where
 * we use a 21-term Taylor series expansion.
 *
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 *
 *      Int_(nu_low)^{nu_high) B(nu,T) dnu
 * 
 * then you must multiply by a factor of acT^4/4pi where a is the radiation
 * constant.  If you want to evaluate expressions like the following:
 *
 *      Int_(4pi)Int_(nu_low)^{nu_high) B(nu,T) dnu domega
 *
 * then you must multiply by acT^4.
 *
 * In the limit of T -> zero, b(T) -> zero, therefore we return a hard zero
 * for a temperature equal to a hard zero.
 *
 * \param lowFreq lower frequency bound in keV
 *
 * \param highFreq higher frequency bound in keV
 *
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Plankian from x_low to x_high
 *
 */
double CDI::integratePlanckSpectrum(double lowFreq, double highFreq, double T)
{
    Require (lowFreq >= 0.0);
    Require (highFreq >= lowFreq);
    Require (T >= 0.0);

    // return 0 if temperature is a hard zero
    if (T == 0.0)
	return 0.0;

    // determine the upper and lower x
    double lower_x = lowFreq  / T;
    double upper_x = highFreq / T;

    // determine the upper and lower bounds calculated by the taylor and
    // polylogarithmic approximations
    double lower_taylor  = taylor_series_planck(lower_x);
    double upper_taylor  = taylor_series_planck(upper_x);
    double lower_poly_m1 = polylog_series_minus_one_planck(lower_x);
    double upper_poly_m1 = polylog_series_minus_one_planck(upper_x);
    double lower_poly    = lower_poly_m1 + 1.0;
    double upper_poly    = upper_poly_m1 + 1.0;

    // determine the integral based on the upper and lower bounds
    double integral = 0.0;

    if (lower_taylor < lower_poly && upper_taylor < upper_poly)
	integral = upper_taylor - lower_taylor;
    else if (lower_taylor < lower_poly && upper_taylor >= upper_poly)
	integral = upper_poly - lower_taylor;
    else 
	integral = upper_poly_m1 - lower_poly_m1;

    Ensure (integral >= 0.0 && integral <= 1.0);

    return integral;
}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the normalized Planckian spectrum from 0 to x (hnu/kT).
 *
 * This function integrates the normalized Plankian that is defined:
 *
 * b(x) = 15/pi^4 * x^3 / (e^x - 1)
 *
 * where 
 *
 * x = (h nu) / (kT)
 * 
 * and 
 * 
 * B(nu, T)dnu = acT^4/4pi * b(x)dx
 * 
 * where B(nu, T) is the Plankian and is defined
 *
 * B(nu, T) = 2hnu^3/c^2 (e^hnu/kt - 1)^-1
 *
 * This function performs the following integration:
 *
 *      Int_(0)^(x) b(x) dx
 *
 * using the method of B. Clark (JCP (70)/2, 1987).  We use a 10-term
 * Polylogarithmic expansion for the normalized Planckian, except in the
 * low-x limit, where we use a 21-term Taylor series expansion.
 *
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 *
 *      Int_(0)^{nu) B(nu,T) dnu
 * 
 * then you must multiply by a factor of acT^4/4pi where a is the radiation
 * constant.  If you want to evaluate expressions like the following:
 *
 *      Int_(4pi)Int_(0)^{nu) B(nu,T) dnu domega
 *
 * then you must multiply by acT^4.
 *
 * In the limit of T -> zero, b(T) -> zero, therefore we return a hard zero
 * for a temperature equal to a hard zero.
 *
 * \param frequency frequency upper integration limit in keV
 *
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Plankian from 0 to x (hnu/kT)
 *
 */
double CDI::integratePlanckSpectrum(double frequency, double T)
{
    Require (T >= 0.0);
    Require (frequency >= 0.0);

    // calculate the integral
    double integral = integratePlanckSpectrum(0.0, frequency, T);
      
    Ensure (integral >= 0.0 && integral <= 1.0);

    return integral;
}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Planckian spectrum over a frequency group.
 *
 * This function integrates the normalized Plankian that is defined:
 *
 * b(x) = 15/pi^4 * x^3 / (e^x - 1)
 *
 * where 
 *
 * x = (h nu) / (kT)
 * 
 * and 
 * 
 * B(nu, T)dnu = acT^4/4pi * b(x)dx
 * 
 * where B(nu, T) is the Plankian and is defined
 *
 * B(nu, T) = 2hnu^3/c^2 (e^hnu/kt - 1)^-1
 *
 * The normalized Plankian, integrated from 0 to infinity, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Planck function may integrate to something less than one.
 *
 * This function performs the following integration:
 *
 *      Int_(x_g-1)^(x_g) b(x) dx
 *
 * using the method of B. Clark (JCP (70)/2, 1987).  We use a 10-term
 * Polylogarithmic expansion for the normalized Planckian, except in the
 * low-x limit, where we use a 21-term Taylor series expansion.
 *
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 *
 *      Int_(nu_g-1)^{nu_g) B(nu,T) dnu
 * 
 * then you must multiply by a factor of acT^4/4pi where a is the radiation
 * constant.  If you want to evaluate expressions like the following:
 *
 *      Int_(4pi)Int_(nu_g-1)^{nu_g) B(nu,T) dnu domega
 *
 * then you must multiply by acT^4.
 *
 * In the limit of T -> zero, b(T) -> zero, therefore we return a hard zero
 * for a temperature equal to a hard zero.
 *
 * If no groups are defined then an exception is thrown.
 *
 * \param groupIndex index of the frequency group to integrate [1,num_groups]
 *
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Plankian over the group specified by
 * groupIndex
 *
 */
double CDI::integratePlanckSpectrum(int groupIndex, double T)
{
    Insist  (!frequencyGroupBoundaries.empty(), "No groups defined!");
    Require (T >= 0.0);
    Require (groupIndex > 0 && 
	     groupIndex <= frequencyGroupBoundaries.size() - 1);

    // first determine the group boundaries for groupIndex
    double lower_bound = frequencyGroupBoundaries[groupIndex-1];
    double upper_bound = frequencyGroupBoundaries[groupIndex];
    Check (upper_bound > lower_bound);

    // calculate the integral over the frequency group
    double integral = integratePlanckSpectrum(lower_bound, upper_bound, T);
	  
    Ensure (integral >= 0.0 && integral <= 1.0);

    return integral;
}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Planckian spectrum over all frequency groups.
 *
 * This function integrates the normalized Plankian that is defined:
 *
 * b(x) = 15/pi^4 * x^3 / (e^x - 1)
 *
 * where 
 *
 * x = (h nu) / (kT)
 * 
 * and 
 * 
 * B(nu, T)dnu = acT^4/4pi * b(x)dx
 * 
 * where B(nu, T) is the Plankian and is defined
 *
 * B(nu, T) = 2hnu^3/c^2 (e^hnu/kt - 1)^-1
 *
 * The normalized Plankian, integrated from 0 to infinity, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Planck function may integrate to something less than one.
 *
 * This function performs the following integration:
 *
 *      Int_(x_1)^(x_N) b(x) dx
 *
 * where x_1 is the low frequency bound and x_N is the high frequency bound
 * of the multigroup data set.  This integration uses the method of B. Clark
 * (JCP (70)/2, 1987).  We use a 10-term Polylogarithmic expansion for the
 * normalized Planckian, except in the low-x limit, where we use a 21-term
 * Taylor series expansion.
 *
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 *
 *      Int_(nu_1)^{nu_N) B(nu,T) dnu
 * 
 * then you must multiply by a factor of acT^4/4pi where a is the radiation
 * constant.  If you want to evaluate expressions like the following:
 *
 *      Int_(4pi)Int_(nu_1)^{nu_N) B(nu,T) dnu domega
 *
 * then you must multiply by acT^4.
 *
 * In the limit of T -> zero, b(T) -> zero, therefore we return a hard zero
 * for a temperature equal to a hard zero.
 *
 * If no groups are defined then an exception is thrown.
 *
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Plankian over all frequency groups
 *
 */
double CDI::integratePlanckSpectrum(double T)
{
    Insist  (!frequencyGroupBoundaries.empty(), "No groups defined!");
    Require (T >= 0.0);

    // first determine the group boundaries for groupIndex
    double lower_bound = frequencyGroupBoundaries.front();
    double upper_bound = frequencyGroupBoundaries.back();
    Check (upper_bound > lower_bound);

    // calculate the integral 
    double integral = integratePlanckSpectrum(lower_bound, upper_bound, T);
	  
    Ensure (integral >= 0.0 && integral <= 1.0);

    return integral;
}

//---------------------------------------------------------------------------//
// SET FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Register a gray opacity (rtt_cdi::GrayOpacity) with CDI.
 *
 * This function sets a gray opacity object of type rtt_cdi::GrayOpacity with
 * the CDI object.  It stores the gray opacity object based upon its
 * rtt_cdi::Model and rtt_cdi::Reaction types.  If a GrayOpacity with these
 * types has already been registered an exception is thrown.  To register a
 * new set of GrayOpacity objects call CDI::reset() first.  You cannot
 * overwrite registered objects with the setGrayOpacity() function!
 *
 * \param spGOp smart pointer (rtt_dsxx::SP) to a GrayOpacity object
 */
void CDI::setGrayOpacity(const SP_GrayOpacity &spGOp)
{
    Require (spGOp);

    // determine the model and reaction type (these MUST be in the correct
    // range because the Model and Reaction are constrained by the
    // rtt_cdi::Model and rtt_cdi::Reaction enumerations, assuming nobody
    // hosed these)
    int model    = spGOp->getModelType();
    int reaction = spGOp->getReactionType();

    Insist (!grayOpacities[model][reaction], 
	    "Tried to overwrite a set GrayOpacity object!");

    // assign the smart pointer
    grayOpacities[model][reaction] = spGOp;

    Ensure (grayOpacities[model][reaction]);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Register a multigroup opacity (rtt_cdi::MultigroupOpacity) with
 * CDI.
 *
 * This function sets a multigroup opacity object of type
 * rtt_cdi::MultigroupOpacity with the CDI object.  It stores the multigroup
 * opacity object based upon its rtt_cdi::Model and rtt_cdi::Reaction types.
 * If a MultigroupOpacity with these type has already been registered an
 * exception is thrown.  To register a new set of MultigroupOpacity objects
 * call CDI::reset() first.  You cannot overwrite registered objects with the
 * setMultigroupOpacity() function!
 *
 * \param spGOp smart pointer (rtt_dsxx::SP) to a MultigroupOpacity object
 */
void CDI::setMultigroupOpacity(const SP_MultigroupOpacity &spMGOp) 
{
    using rtt_dsxx::soft_equiv;

    Require (spMGOp);

    // determine the model and reaction types
    int model    = spMGOp->getModelType();
    int reaction = spMGOp->getReactionType();

    Insist (!multigroupOpacities[model][reaction],
	    "Tried to overwrite a set MultigroupOpacity object!");

    // if the frequency group boundaries have not been assigned in any CDI
    // object, then assign them here
    if (frequencyGroupBoundaries.empty())
    {
	// copy the the group boundaries for this material to the "global"
	// group boundaries that will be enforced for all CDI objects
	frequencyGroupBoundaries = spMGOp->getGroupBoundaries();
    }

    // always check that the number of frequency groups is the same for each
    // multigroup material added to CDI
    Insist (spMGOp->getNumGroupBoundaries() == 
	    frequencyGroupBoundaries.size(),
	    "Incompatible frequency groups assigned for this material");

    // do a check of the actual boundary values when DBC check is on (this is
    // more expensive so we retain the option of turning it off)
    const std::vector<double> &ref = spMGOp->getGroupBoundaries();
    Check (soft_equiv(frequencyGroupBoundaries.begin(),
		      frequencyGroupBoundaries.end(),
		      ref.begin(),
		      ref.end(),
		      1.0e-6));

    // assign the smart pointer
    multigroupOpacities[model][reaction] = spMGOp;

    Ensure (multigroupOpacities[model][reaction]);
}

//---------------------------------------------------------------------------//

void CDI::setEoS( const SP_EoS &in_spEoS )
{
    Require (in_spEoS);

    Insist (!spEoS, "Tried to overwrite a set EoS object.!");

    // set the smart pointer
    spEoS = in_spEoS;

    Ensure (spEoS);
}

//---------------------------------------------------------------------------//
// GET FUNCTIONS
//---------------------------------------------------------------------------//

// Provide CDI with access to the full interfaces defined by GrayOpacity.hh
// and MultigroupOpacity.hh
    
CDI::SP_GrayOpacity CDI::gray(rtt_cdi::Model    m, 
			      rtt_cdi::Reaction r) const 
{ 
    Insist (grayOpacities[m][r], "Undefined GrayOpacity!");
    return grayOpacities[m][r]; 
}
    
CDI::SP_MultigroupOpacity CDI::mg(rtt_cdi::Model    m,
				  rtt_cdi::Reaction r) const 
{
    Insist (multigroupOpacities[m][r], "Undefined MultigroupOpacity!");
    return multigroupOpacities[m][r]; 
}

//---------------------------------------------------------------------------//

// Provide CDI with access to the full interfaces defined by
// EoS.hh
    
CDI::SP_EoS CDI::eos() const 
{ 
    Insist (spEoS, "Undefined EoS!");
    return spEoS; 
}

//---------------------------------------------------------------------------//
// RESET THE CDI OBJECT
//---------------------------------------------------------------------------//

void CDI::reset()
{
    Check (grayOpacities.size() == constants::num_Models);
    Check (multigroupOpacities.size() == constants::num_Models);

    // reset the gray opacities
    for (int i = 0; i < constants::num_Models; i++)
    {
	Check (grayOpacities[i].size() == constants::num_Reactions);
	Check (multigroupOpacities[i].size() == constants::num_Reactions);

	for (int j = 0; j < constants::num_Reactions; j++)
	{
	    // reassign the GrayOpacity SP to a null SP
	    grayOpacities[i][j] = SP_GrayOpacity();

	    // reassign the MultigroupOpacity SP to a null SP
	    multigroupOpacities[i][j] = SP_MultigroupOpacity();

	    // check it
	    Check (!grayOpacities[i][j]);
	    Check (!multigroupOpacities[i][j]);
	}
    }

    // empty the frequency group boundaries
    frequencyGroupBoundaries.clear();
    Check (frequencyGroupBoundaries.empty());

    // reset the EoS SP
    spEoS = SP_EoS();
    Check (!spEoS);
}

//---------------------------------------------------------------------------//
// BOOLEAN QUERY FUNCTIONS
//---------------------------------------------------------------------------//

bool CDI::isGrayOpacitySet(rtt_cdi::Model m, rtt_cdi::Reaction r) const
{
    return static_cast<bool>(grayOpacities[m][r]);
}


bool CDI::isMultigroupOpacitySet(rtt_cdi::Model m, rtt_cdi::Reaction r) const
{
    return static_cast<bool>(multigroupOpacities[m][r]);
}

bool CDI::isEoSSet() const
{
    return static_cast<bool>(spEoS);
}

} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
//                              end of CDI.cc
//---------------------------------------------------------------------------//
