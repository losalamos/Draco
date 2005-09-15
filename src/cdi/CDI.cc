//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:07 2000
 * \brief  CDI class implementation file.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "CDI.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"
#include <cmath>
#include <limits>

namespace rtt_cdi
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS AND DESTRUCTORS
//---------------------------------------------------------------------------//

/*!
 * \brief Construct a CDI object.
 *
 * Builds a CDI object.  The opacity and eos objects that this holds must
 * be loaded using the set functions.  There is no easy way to guarantee
 * that all of the set objects point to the same material.  CDI does do
 * checking that only one of each Model:Reaction pair of opacity objects
 * are assigned; however, the user can "fake" CDI with different
 * materials if he/she is malicious enough.  
 *
 * CDI does allow a string material ID indicator.  It is up to the client
 * to ascribe meaning to the indicator.
 *
 * \param id string material id descriptor, this is defaulted to null
 */
CDI::CDI(const std_string &id)
    : matID(id),
      grayOpacities(
          constants::num_Models, 
          SF_GrayOpacity(constants::num_Reactions)),
      multigroupOpacities(
          constants::num_Models,
          SF_MultigroupOpacity(constants::num_Reactions))
{

    Ensure (grayOpacities.size() == constants::num_Models);
    Ensure (multigroupOpacities.size() == constants::num_Models);

}

//---------------------------------------------------------------------------//
    
CDI::~CDI() { /* empty */ }


//---------------------------------------------------------------------------//
// STATIC DATA
//---------------------------------------------------------------------------//

std::vector<double> CDI::frequencyGroupBoundaries = std::vector<double>();

//---------------------------------------------------------------------------//
// STATIC FUNCTIONS
//---------------------------------------------------------------------------//

/*!
 * \brief Return the frequency group boundaries.
 *
 * Every multigroup opacity object held by any CDI object contains the
 * same frequency group boundaries.  This static function allows CDI
 * users to access the group boundaries without referencing a particular
 * material. 
 *
 * Note, the group boundaries are not set until a multigroup opacity
 * object is set for the first time (in any CDI object) with the
 * setMultigroupOpacity function.
 */
std::vector<double> CDI::getFrequencyGroupBoundaries()
{
    return frequencyGroupBoundaries;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the number of frequency groups.
 */
int CDI::getNumberFrequencyGroups()
{
    int ng = 0;

    if (!frequencyGroupBoundaries.empty())
	ng = frequencyGroupBoundaries.size() - 1;

    return ng;
}    

//---------------------------------------------------------------------------//
// UNNAMED NAMESPACE
//---------------------------------------------------------------------------//
// Nested unnamed namespace that holds data and services used by the
// Planckian integration routines.  The data in this namespace is accessible
// by the methods in this file only (internal linkage).

namespace
{

// Constants used in the Taylor series expansion of the Planckian:
static const double coeff_3  =    1.0 / 3.0;
static const double coeff_4  =   -1.0 / 8.0;
static const double coeff_5  =    1.0 / 60.0;
static const double coeff_7  =   -1.0 / 5040.0;
static const double coeff_9  =    1.0 / 272160.0;
static const double coeff_11 =   -1.0 / 13305600.0;
static const double coeff_13 =    1.0 / 622702080.0;
static const double coeff_15 =  -6.91 / 196151155200.0;
static const double coeff_17 =    1.0 / 1270312243200.0;
static const double coeff_19 = -3.617 / 202741834014720.0;
static const double coeff_21 = 43.867 / 107290978560589824.0;

static const double coeff    =   0.1539897338202651; // 15/pi^4

//---------------------------------------------------------------------------//
/*!
 * \brief Computes the normalized Planck integral via a 21 term Taylor
 * expansion. 
 *
 * The taylor expansion of the planckian integral looks as follows:
 *
 * I(x) = c0(c3x^3 + c4x^4 + c5x^5 + c7x^7 + c9x^9 + c11x^11 + c13x^13 +
 *           c15x^15 + c17x^17 + c19x^19 + c21x^21)
 *
 * If done naively, this requires 136 multiplications. If you accumulate
 * the powers of x as you go, it can be done with 24 multiplications.
 *
 * If you express the polynomaial as follows:
 * 
 * I(x) = c0*x^3(c3 + x(c4 + x(c5 + x^2(c7 + x^2(c9 + x^2(c11 + x^2(c13 +
 *               x^2(c15 + x^2(c17 + x^2(c19 + x^2c21))))))))))
 *
 * the evaluation can be done with 13 multiplications. Furthermore, we do
 * not need to worry about overflow on large powers of x, since the largest
 * power we compute is x^3
 *
 * \param  The point at which the Planck integral is evaluated.
 * \return The integral value.
 */

inline double taylor_series_planck(double x)
{

    Require( x >= 0.0 );

    const double xsqrd  = x * x;

    double taylor ( coeff_21 * xsqrd );

    taylor += coeff_19;
    taylor *= xsqrd;

    taylor += coeff_17;
    taylor *= xsqrd;

    taylor += coeff_15;
    taylor *= xsqrd;

    taylor += coeff_13;
    taylor *= xsqrd;

    taylor += coeff_11;
    taylor *= xsqrd;

    taylor += coeff_9;
    taylor *= xsqrd;

    taylor += coeff_7;
    taylor *= xsqrd;

    taylor += coeff_5;
    taylor *= x;

    taylor += coeff_4;
    taylor *= x;

    taylor += coeff_3;
    taylor *= x * xsqrd * coeff;
    
    Ensure (taylor >= 0.0);

    return taylor;
}


// ---------------------------------------------------------------------------
// return the 10-term Polylogarithmic expansion (minus one) for the Planck
// integral given x
double polylog_series_minus_one_planck(const double x)
{
    Require (x >= 0.0);

    const double eix = std::exp(-x);
    double xsqrd = x * x;

    static const double i_plus_two_inv[9] = 
	{
	    0.5000000000000000, // 1/2
	    0.3333333333333333, // 1/3
	    0.2500000000000000, // 1/4
	    0.2000000000000000, // 1/5
	    0.1666666666666667, // 1/6
	    0.1428571428571429, // 1/7
	    0.1250000000000000, // 1/8
	    0.1111111111111111, // 1/9
	    0.1000000000000000  // 1/10
	};
    double const * curr_inv = i_plus_two_inv;
    
    // initialize to what would have been the calculation of the i=1 term.
    // This saves a number of "mul by one" ops.
    double eixp = eix;
    double li1  = eix;
    double li2  = eix; 
    double li3  = eix;
    double li4  = eix;

    // calculate terms 2..10.  This loop has been unrolled by a factor of 3
    for(int i = 2; i < 11; i += 3)
    {
	register const double ip0_inv = *curr_inv++;
	eixp *= eix;
	double eixr_ip0 = eixp * ip0_inv;

	register const double ip1_inv = *curr_inv++;
	eixp *= eix;
	double eixr_ip1 = eixp * ip1_inv;

	register const double ip2_inv = *curr_inv++;
	eixp *= eix;
	double eixr_ip2 = eixp * ip2_inv;

	
	const double r10 = eixr_ip0;
	const double r11 = eixr_ip1;
	const double r12 = eixr_ip2;

	const double r20 = (eixr_ip0 *= ip0_inv);
	const double r21 = (eixr_ip1 *= ip1_inv);
	const double r22 = (eixr_ip2 *= ip2_inv);

	li1 += r10 + r11 + r12;

	const double r30 = (eixr_ip0 *= ip0_inv);
	const double r31 = (eixr_ip1 *= ip1_inv);
	const double r32 = (eixr_ip2 *= ip2_inv);

	li2 += r20 + r21 + r22;

	const double r40 = (eixr_ip0 *= ip0_inv);
	const double r41 = (eixr_ip1 *= ip1_inv);
	const double r42 = (eixr_ip2 *= ip2_inv);

	li3 += r30 + r31 + r32;
	li4 += r40 + r41 + r42;
    }



    register const double xcubed = xsqrd*x;
    xsqrd *= 3.0;			// <<<<<<< LOOK: now 3x^2!!!!!!!!

    // calculate the lower polylogarithmic integral
    const double poly = -coeff * (xcubed * li1 +
				  xsqrd  * li2 + // really, xsqrd*3
				  6.0 * (x * li3 + li4));
    
    Ensure (poly <= 0.0);
    return poly;
}



double Planck2Rosseland(const double freq)
{

    static const double NORM_FACTOR = 0.25*coeff; //  15./(4.*PI4);

    double e_freq = exp(-freq);
    double freq_3 = freq*freq*freq;
    
    double factor;

    if (freq > 1.0e-5)
        factor = NORM_FACTOR * e_freq * (freq_3*freq) / (1 - e_freq);
    else
        factor = NORM_FACTOR * freq_3 / (1 - 0.5*freq);

    return factor;

}

} // end of unnamed namespace


//---------------------------------------------------------------------------//
// Core Integrators
/*
 * These are the most basic of the Planckian and Rosseland integration
 * functions. They are publically accessible, but also used in the
 * implementation of integration functions with friendlier interfaces.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Integrate the normalized Planckian spectrum from 0 to \f$ x
 * (\frac{h\nu}{kT}) \f$.
 *
 * \param freq frequency upper integration limit in keV
 *
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Plankian from 0 to x \f$(\frac{h\nu}{kT})\f$
 *
 *
 */
double CDI::integratePlanckSpectrum(const double freq, 
                                    const double T) 

{
    Require (freq >= 0);
    Require (T    >= 0.0);

    // Return 0 if temperature is a hard zero
    if (T == 0.0) return 0.0;

    const double scaled = freq / T;
    const double poly   = polylog_series_minus_one_planck(scaled) + 1.0;
    const double taylor = taylor_series_planck(scaled);

    double integral = std::min(taylor, poly);

    Ensure ( integral >= 0.0 );
    Ensure ( integral <= 1.0 );

    return integral;
}


//---------------------------------------------------------------------------//
/* \brief
 *
 */
double CDI::integrateRosselandSpectrum(const double freq,
                                       const double T)
{

    Require (freq >= 0.0);
    Require (T    >= 0.0);

    double planck, rosseland;

    integratePlanckRosselandSpectrum(freq, T, planck, rosseland);

    return rosseland;

}


//---------------------------------------------------------------------------//
/* \brief
 *
 */
void CDI::integratePlanckRosselandSpectrum(const double freq,
                                           const double T,
                                           double& planck,
                                           double& rosseland)
{

    Require (freq >= 0.0);
    Require (T    >= 0.0);

    // Return 0 if temperature is a hard zero.
    if (T == 0.0)
    {
        planck    = 0;
        rosseland = 0;
	return;
    }

    // Calculate the Planckian integral 
    planck = integratePlanckSpectrum(freq, T);
    
    Ensure (planck >= 0.0);
    Ensure (planck <= 1.0);

    double scaled        = freq / T;
    double add_rosseland = Planck2Rosseland(scaled);

    rosseland = planck - add_rosseland;

}



//---------------------------------------------------------------------------//
// Planckian Spectrum Integrators
//
/* These are versions of the integrators that work over specific energy ranges
 * or groups in the stored group structure.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Planckian spectrum over a frequency range.
 *
 * \param lowFreq lower frequency bound in keV
 * \param highFreq higher frequency bound in keV
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Plankian from x_low to x_high
 *
 */
double CDI::integratePlanckSpectrum(const double lowFreq, 
				    const double highFreq, 
				    const double T) 
{
    Require (lowFreq  >= 0.0);
    Require (highFreq >= lowFreq);
    Require (T        >= 0.0);

    // return 0 if temperature is a hard zero
    if (T == 0.0) return 0.0;

    double integral =
        integratePlanckSpectrum(highFreq,T) -
        integratePlanckSpectrum(lowFreq,T);
    

    Ensure ( integral >= 0.0 );
    Ensure ( integral <= 1.0 );

    return integral;
}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Planckian spectrum over a frequency group.
 *
 * \param groupIndex Index of the frequency group to integrate [1,num_groups].
 * \param T          The temperature in keV (must be greater than 0.0).
 * \return           Integrated normalized Plankian over the group specified
 *                   by groupIndex.
 *
 */
double CDI::integratePlanckSpectrum(const int groupIndex, const double T)
{
    Insist  (!frequencyGroupBoundaries.empty(), "No groups defined!");

    Require (T >= 0.0);
    Require (groupIndex > 0);
    Require (groupIndex <= frequencyGroupBoundaries.size() - 1);

    // Determine the group boundaries for groupIndex
    double lower_bound = frequencyGroupBoundaries[groupIndex-1];
    double upper_bound = frequencyGroupBoundaries[groupIndex];
    Check (upper_bound > lower_bound);

    double integral = integratePlanckSpectrum(lower_bound, upper_bound, T);
	  
    Ensure (integral >= 0.0);
    Ensure (integral <= 1.0);

    return integral;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Integrate the Planckian spectrum over all frequency groups.
 * \param T The temperature in keV (must be greater than 0.0).
 * \return Integrated normalized Plankian over all frequency groups.
 */
double CDI::integratePlanckSpectrum(const double T)
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
// Rosseland Spectrum Integrators
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Integrate the Planckian and Rosseland spectrum over a frequency
 * range.
 * \param T the temperature in keV (must be greater than 0.0)
 * \return void the integrated normalized Planckian and Rosseland from x_low
 * to x_high are passed by reference in the function call
 *
 */
void CDI::integrate_Rosseland_Planckian_Spectrum(const double  lowFreq,
						 const double  highFreq,
						 const double  T,
						 double       &planck, 
						 double       &rosseland)
{
    Require (lowFreq >= 0.0);
    Require (highFreq >= lowFreq);
    Require (T >= 0.0);

    double planck_high, rosseland_high;
    double planck_low,  rosseland_low;

    integratePlanckRosselandSpectrum(lowFreq,  T, planck_low,  rosseland_low);
    integratePlanckRosselandSpectrum(highFreq, T, planck_high, rosseland_high);

    planck    = planck_high    - planck_low;
    rosseland = rosseland_high - rosseland_low;

}
    


//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Rosseland spectrum over a frequency range.
 *
 
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Rosseland from x_low to x_high
 *
 */
double CDI::integrateRosselandSpectrum(const double lowFreq,
				       const double highFreq, 
				       const double T)
{

    Require (lowFreq  >= 0.0);
    Require (highFreq >= lowFreq);
    Require (T >= 0.0);

    double planck, rosseland;

    integrate_Rosseland_Planckian_Spectrum(lowFreq, highFreq, T, planck, rosseland); 

    return rosseland;

}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Rosseland spectrum over a frequency group.
 *
 * \param groupIndex index of the frequency group to integrate [1,num_groups]
 * \param T          the temperature in keV (must be greater than 0.0)
 * \return           integrated normalized Plankian over the group specified
 *                   by groupIndex.
 *
 */
double CDI::integrateRosselandSpectrum(const int groupIndex, const double T)
{
    Insist  (!frequencyGroupBoundaries.empty(), "No groups defined!");
    Require (T >= 0.0);
    Require (groupIndex > 0 && 
	     groupIndex <= frequencyGroupBoundaries.size() - 1);

    // first determine the group boundaries for groupIndex
    double lowFreq  = frequencyGroupBoundaries[groupIndex-1];
    double highFreq = frequencyGroupBoundaries[groupIndex];
    Check (highFreq  > lowFreq);

    double rosseland = integrateRosselandSpectrum(lowFreq, highFreq, T);

    return rosseland;

}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Planckian and Rosseland spectrum over a frequency
 * group.
 *
 * \param groupIndex index of the frequency group to integrate [1,num_groups]
 * \param T          The temperature in keV (must be greater than 0.0)
 * \param PL         Reference argument for the Planckian integral
 * \param ROSL       Reference argument for the Rosseland integral
 *
 * \return The integrated normalized Planckian and Rosseland over the
 * requested frequency group. These are returned as references in argument PL
 * and ROSL
 *
 */
void CDI::integrate_Rosseland_Planckian_Spectrum(const int     groupIndex,
						 const double  T, 
						 double       &planck, 
						 double       &rosseland)
{
    Insist  (!frequencyGroupBoundaries.empty(), "No groups defined!");

    Require (T >= 0.0);
    Require (groupIndex > 0);
    Require (groupIndex <= frequencyGroupBoundaries.size() - 1);

    // Determine the group boundaries
    double lowFreq  = frequencyGroupBoundaries[groupIndex-1];
    double highFreq = frequencyGroupBoundaries[groupIndex];
    Check (highFreq  > lowFreq);

    // Call the general frequency version
    integrate_Rosseland_Planckian_Spectrum(lowFreq, highFreq, T, planck, rosseland);
 
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

/*!
 * \brief This fuction returns a GrayOpacity object.
 *
 * This provides the CDI with the full functionality of the interface
 * defined in GrayOpacity.hh.  For example, the host code could make the
 * following call: <tt> double newOp = spCDI1->gray()->getOpacity(
 * 55.3, 27.4 ); </tt>
 *
 * The appropriate gray opacity is returned for the given model and
 * reaction type.
 *
 * \param m rtt_cdi::Model specifying the desired physics model
 * \param r rtt_cdi::Reaction specifying the desired reaction type
 */
CDI::SP_GrayOpacity CDI::gray(rtt_cdi::Model    m, 
			      rtt_cdi::Reaction r) const 
{ 
    Insist (grayOpacities[m][r], "Undefined GrayOpacity!");
    return grayOpacities[m][r]; 
}
    
/*!
 * \brief This fuction returns the MultigroupOpacity object.
 *
 * This provides the CDI with the full functionality of the interface
 * defined in MultigroupOpacity.hh.  For example, the host code could
 * make the following call:<br> <tt> int numGroups =
 * spCDI1->mg()->getNumGroupBoundaries(); </tt>
 *
 * The appropriate multigroup opacity is returned for the given reaction
 * type.
 *
 * \param m rtt_cdi::Model specifying the desired physics model
 * \param r rtt_cdi::Reaction specifying the desired reaction type.
 *
 */
CDI::SP_MultigroupOpacity CDI::mg(rtt_cdi::Model    m,
				  rtt_cdi::Reaction r) const 
{
    Insist (multigroupOpacities[m][r], "Undefined MultigroupOpacity!");
    return multigroupOpacities[m][r]; 
}

//---------------------------------------------------------------------------//

/*!
 * \brief This fuction returns the EoS object.
 *
 * This provides the CDI with the full functionality of the interface
 * defined in EoS.hh.  For example, the host code could make the
 * following call:<br> <tt> double Cve =
 * spCDI1->eos()->getElectronHeatCapacity( * density, temperature );
 * </tt>
 */
CDI::SP_EoS CDI::eos() const 
{ 
    Insist (spEoS, "Undefined EoS!");
    return spEoS; 
}

//---------------------------------------------------------------------------//
// RESET THE CDI OBJECT
//---------------------------------------------------------------------------//

/*!
 * \brief Reset the CDI object.
 *
 * This function "clears" all data objects (GrayOpacity,
 * MultigroupOpacity, EoS) held by CDI.  After clearing, new objects can
 * be set using the set functions.  
 *
 * As stated in the set functions documentation, you are not allowed to
 * overwrite a data object with the same attributes as one that already
 * has been set.  The only way to "reset" these objects is to call
 * CDI::reset().  Note that CDI::reset() resets \b ALL of the objects
 * stored by CDI (including group boundaries).
 */
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

/*!
 * \brief Query to see if a gray opacity is set.
 */
bool CDI::isGrayOpacitySet(rtt_cdi::Model m, rtt_cdi::Reaction r) const
{
    return static_cast<bool>(grayOpacities[m][r]);
}


/*!
 * \brief Query to see if a multigroup opacity is set.
 */
bool CDI::isMultigroupOpacitySet(rtt_cdi::Model m, rtt_cdi::Reaction r) const
{
    return static_cast<bool>(multigroupOpacities[m][r]);
}

/*!
 * \brief Query to see if an eos is set.
 */
bool CDI::isEoSSet() const
{
    return static_cast<bool>(spEoS);
}

} // end namespace rtt_cdi

//---------------------------------------------------------------------------//
//                              end of CDI.cc
//---------------------------------------------------------------------------//
