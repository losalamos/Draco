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
#include <limits>

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
    using std::numeric_limits;
    using std::log;

    Require( x >= 0.0 );

    // Check for potential overflow errors:

    if( x > 1.0 )
    {
	// This Taylor series expansion takes x to the 21st power.  We need
	// to ensure that x is small enough so that x^21 is less than the
	// largest double that can be represented on the current machine.

	// The maximum double for the current machine is
	double const maxDouble( numeric_limits<double>::max() );

	// To be conservative, assume that we are looking for x^22 so that we
	// will now require that x^22 < maxDouble.  However, if x is too big
	// this comparison will also result in an overflow.  If we take the
	// log of both sides, we will have:
	//
	// Require( log( x^22 ) < log(maxDouble) );
	//
	// This can be simplified to avoid the overflow as follows:
	
	Require( 22.0*log(x) < log(maxDouble) );
    }

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
 * \f[
 *    b(x) = \frac{15}{\pi^4} \frac{x^3}{e^x - 1}
 * \f]
 *
 * where 
 * 
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 *
 * and 
 *
 * \f[ 
 *    B(\nu,T)d\nu = \frac{acT^4}{4\pi} b(x)dx
 * \f]
 * 
 * where \f$B(\nu,T)\f$ is the Plankian and is defined
 *
 * \f[
 *    B(\nu,T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{h\nu/kt} - 1}
 * \f]
 *
 * The normalized Plankian, integrated from 0 to \f$\infty\f$, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Planck function may integrate to something less than one.
 *
 * This function performs the following integration:
 *\f[
 *      \int_{x_low}^{x_high} b(x) dx
 *\f]
 * where \f$x_low\f$ is calculated from the input low frequency bound and
 * \f$x_high\f$ is calculated from the input high frequency bound.  This
 * integration uses the method of B. Clark (JCP (70)/2, 1987).  We use a
 * 10-term Polylogarithmic expansion for the normalized Planckian, except in
 * the low-x limit, where we use a 21-term Taylor series expansion.
 *
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 * \f[
 *      \int_{\nu_low}^{\nu_high} B(\nu,T) d\nu 
 * \f]
 * then you must multiply by a factor of \f$\frac{acT^4}{4\pi}\f$ where a is
 * the radiation constant.  If you want to evaluate expressions like the
 * following:
 *\f[
 *      \int_{4\pi} \int_{\nu_low}^{\nu_high} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$acT^4\f$.
 *
 * In the limit of \f$T \rightarrow 0, b(T) \rightarrow 0, therefore we
 * return a hard zero for a temperature equal to a hard zero.
 *
 * The integral is calculated using a polylogarithmic series approximation,
 * except for low frequencies, where the Taylor series approximation is more
 * accurate.  Each of the Taylor and polylogarithmic approximation has a
 * positive trucation error, so they intersect above the correct solution;
 * therefore, we always use the smaller one for a continuous concatenated
 * function.  When both frequency bounds reside above the Planckian peak
 * (above 2.822 T), we skip the Taylor series calculations and use the
 * polylogarithmic series minus one (the minus one is for roundoff control).
 *
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
double CDI::integratePlanckSpectrum(const double lowFreq, 
				    const double highFreq, 
				    const double T) 
{
    Require (lowFreq  >= 0.0);
    Require (highFreq >= lowFreq);
    Require (T        >= 0.0);

    // return 0 if temperature is a hard zero
    if (T == 0.0)
	return 0.0;

    // determine the upper and lower x
    const double lower_x = lowFreq  / T;
    const double upper_x = highFreq / T;

    // initialize the return integral value
    double integral = 0.0;

    // determine the upper and lower bounds calculated by the 
    // polylogarithmic approximations minus one
    const double lower_poly_m1 = polylog_series_minus_one_planck(lower_x);
    const double upper_poly_m1 = polylog_series_minus_one_planck(upper_x);

    // if both reduced frequencies are above the Planckian peak, 2.82144,
    // bypass the Taylor series approximation and use only the
    // polylogarithmic approximations minus one.
    const double x_at_planck_peak = 2.822;

    if (lower_x >= x_at_planck_peak)
    {
	integral = upper_poly_m1 - lower_poly_m1;
    }

    // otherwise, calculate the taylor and polylogarithmic approximations for
    // the upper and lower bounds and use the appropriate ones.  Both have
    // positive truncation errors, and they intersect above the true
    // function, therefore we always use the minimum.  Using the minimum
    // means that we do not have to explicitly calculate the intersection.
    else
    {
	const double lower_taylor = taylor_series_planck(lower_x);
	const double upper_taylor = taylor_series_planck(upper_x);
	const double lower_poly   = lower_poly_m1 + 1.0;
	const double upper_poly   = upper_poly_m1 + 1.0;

	// both limits are below the intersection, so use Taylor for both
	if ( upper_taylor < upper_poly )
	{
	    Check ( lower_taylor < lower_poly );
	    integral = upper_taylor - lower_taylor;
	}

	// the limits straddle the intersection, use both Taylor and polylog
	else if ( lower_taylor < lower_poly )
	{
	    Check ( upper_taylor >= upper_poly );
	    integral = upper_poly - lower_taylor;
	}

	// both limits are above the intersection, so use polylog-1 for both
	else 
	{
	    Check ( lower_taylor >= lower_poly );
	    integral = upper_poly_m1 - lower_poly_m1;
	}
    }

    Ensure ( integral >= 0.0 )
    Ensure ( integral <= 1.0 );

    return integral;
}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the normalized Planckian spectrum from 0 to \f$ x
 * (\frac{h\nu}{kT}) \f$.
 *
 * This function integrates the normalized Plankian that is defined:
 *\f[
 *    b(x) = \frac{15}{\pi^4} \frac{x^3}{e^x - 1}
 *\f]
 * where 
 *\f[
 * x = \frac{h\nu}{kT}
 * \f]
 * and 
 * \f[
 *    B(\nu,T)d\nu = \frac{acT^4}{4\pi} b(x)dx
 * \f]
 * where \f$B(\nu,T)\f$ is the Plankian and is defined
 *\f[
 *    B(\nu,T) = \frac{2hnu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}} - 1}
 *\f]
 * This function performs the following integration:
 *\f[
 *      \int_{0}^{x} b(x) dx
 *\f]
 * using the method of B. Clark (JCP (70)/2, 1987).  We use a 10-term
 * Polylogarithmic expansion for the normalized Planckian, except in the
 * low-x limit, where we use a 21-term Taylor series expansion.
 *
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 *\f[
 *     \int_{0}^{\nu} B(\nu,T) d\nu
 * \f]
 * then you must multiply by a factor of \f$\frac{acT^4}{4\pi}\f$ where a is
 * the radiation constant.  If you want to evaluate expressions like the
 * following: 
 *\f[
 *     \int_{4\pi} \int_{0}^{\nu} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$acT^4\f$.
 *
 * In the limit of \f$T \rightarrow 0, b(T) \rightarrow 0 \f$, therefore we
 * return a hard zero for a temperature equal to a hard zero.
 *
 * \param frequency frequency upper integration limit in keV
 *
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Plankian from 0 to x \f$(\frac{h\nu}{kT})\f$
 *
 */
double CDI::integratePlanckSpectrum(const double frequency, const double T)
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
 * \f[
 *    b(x) = \frac{15}{\pi^4} \frac{x^3}{e^x - 1}
 * \f]
 * where 
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 * and 
 * \f[
 *    B(\nu, T)d\nu = \frac{acT^4}{4\pi}b(x)dx
 * \f]
 * where \f$B(\nu, T)\f$ is the Plankian and is defined
 * \f[
 *    B(\nu, T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}} - 1}
 * \f]
 * The normalized Plankian, integrated from 0 to \f$\infty\f$, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Planck function may integrate to something less than one.
 *
 * This function performs the following integration:
 * \f[
 *      \int_{x_{g-1}}^{x_g} b(x) dx
 * \f]
 * using the method of B. Clark (JCP (70)/2, 1987).  We use a 10-term
 * Polylogarithmic expansion for the normalized Planckian, except in the
 * low-x limit, where we use a 21-term Taylor series expansion.
 *
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 * \f[
 *      \int_{\nu_{g-1}}^{\nu_g} B(\nu,T) d\nu
 * \f]
 * then you must multiply by a factor of \f$\frac{acT^4}{4\pi}\f$ where a is
 * the radiation constant.  If you want to evaluate expressions like the
 * following:
 *\f[
 *     \int_{4\pi} \int_{\nu_{g-1}}^{\nu_g} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$acT^4\f$.
 *
 * In the limit of \f$T \rightarrow 0, b(T) \rightarrow 0\f$, therefore we
 * return a hard zero for a temperature equal to a hard zero.
 *
 * If no groups are defined then an exception is thrown.
 *
 * \param groupIndex index of the frequency group to integrate [1,num_groups]
 * \param T          the temperature in keV (must be greater than 0.0)
 * \return           integrated normalized Plankian over the group specified
 *                   by groupIndex.
 *
 */
double CDI::integratePlanckSpectrum(const int groupIndex, const double T)
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
 * \f[
 *    b(x) = \frac{15}{\pi^4} \frac{x^3}{e^x - 1}
 * \f]
 * where 
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 * and 
 * \f[
 *    B(\nu, T)d\nu = \frac{acT^4}{4\pi}b(x)dx
 * \f]
 * where \f$B(\nu, T)\f$ is the Plankian and is defined
 * \f[
 *    B(\nu, T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}} - 1}
 * \f]
 * The normalized Plankian, integrated from 0 to \f$\infty\f$, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Planck function may integrate to something less than one.
 *
 * This function performs the following integration:
 * \f[
 *      \int_{x_1}^{x_N} b(x) dx
 * \f]
 * where \f$x_1\f$ is the low frequency bound and \f$x_N\f$ is the high
 * frequency bound of the multigroup data set.  This integration uses the
 * method of B. Clark (JCP (70)/2, 1987).  We use a 10-term Polylogarithmic
 * expansion for the normalized Planckian, except in the low-x limit, where
 * we use a 21-term Taylor series expansion.
 *
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 * \f[
 *      \int_{\nu_1}^{\nu_N} B(\nu,T) d\nu
 * \f]
 * then you must multiply by a factor of \f$\frac{acT^4}{4\pi}\f$ where a is
 * the radiation constant.  If you want to evaluate expressions like the
 * following:
 *\f[
 *      \int_{4\pi} \int_{\nu_1}^{\nu_N} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$ acT^4 \f$.
 *
 * In the limit of \f$ T \rightarrow 0, b(T) \rightarrow 0\f$, therefore we
 * return a hard zero for a temperature equal to a hard zero.
 *
 * If no groups are defined then an exception is thrown.
 *
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Plankian over all frequency groups
 *
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
/*!
 *
 * \brief Integrate the Rosseland spectrum over a frequency group.
 *
 * This function integrates the normalized Rosseland that is defined:
 * \f[
 *    r(x) = \frac{15}{4\pi^4} \frac{x^4 e^x}{(e^x - 1)^2}
 * \f]
 * where 
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 * and 
 * \f[
 *    R(\nu, T)d\nu = \frac{4 acT^3}{4\pi}r(x)dx
 * \f]
 * where \f$R(\nu, T)\f$ is the Rosseland and is defined
 * \f[
 *    R(\nu, T) = \frac{\partial B(\nu, T)}{\partial T}
 * \f]
 * \f[
 *    B(\nu, T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}} - 1}
 * \f]
 * The normalized Rosseland, integrated from 0 to \f$\infty\f$, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Rosseland function may integrate to something less than one.
 *
 * This function performs the following integration:
 * \f[
 *      \int_{x_1}^{x_N} r(x) dx
 * \f]
 * where \f$x_1\f$ is the low frequency bound and \f$x_N\f$ is the high
 * frequency bound of the multigroup data set.  This integration uses the
 * method of B. Clark (JCP (70)/2, 1987).  We use a 10-term Polylogarithmic
 * expansion for the normalized Planckian, except in the low-x limit, where
 * we use a 21-term Taylor series expansion.
 *
 * For the Rosseland we can relate the group interval integration to the
 * Planckian group interval integration, by equation 27 in B. Clark paper.
 * \f[
 *     \int_{x_1}^{x_N} r(x) dx = \int_{x_1}^{x_N} b(x) dx
 *     - \frac{15}{4\pi^4} \frac{x^4}{e^x - 1}
 * \f]
 * Therefore our Rosslenad group integration function can simply wrap the
 * Planckian function integratePlanckSpectrum(lowFreq, highFreq, T)
 * 
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 * \f[
 *      \int_{\nu_1}^{\nu_N} R(\nu,T) d\nu
 * \f]
 * then you must multiply by a factor of \f$\frac{4 acT^3}{4\pi}\f$ where a is
 * the radiation constant.  If you want to evaluate expressions like the
 * following:
 *\f[
 *      \int_{4\pi} \int_{\nu_1}^{\nu_N} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$ 4 acT^3 \f$.
 *
 * In the limit of \f$ T \rightarrow 0, r(T) \rightarrow 0\f$, therefore we
 * return a hard zero for a temperature equal to a hard zero.
 
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Rosseland from x_low to x_high
 *
 */

double CDI::integrateRosselandSpectrum(const double lowFreq,const double highFreq, const double T)
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
    double NORM_FACTOR = 0.25*coeff; //  15./(4.*PI4);

    // calculate the Planckian integral 
    double integral = integratePlanckSpectrum(lowFreq, highFreq, T);
    
    // Calculate the addition to the Planck integral 
    // for x > xlim=1.e-5
    //  double add_lower = NORM_FACTOR*pow(lower_x,4)/(exp(lower_x)-1.);
    //  double add_upper = NORM_FACTOR*pow(upper_x,4)/(exp(upper_x)-1.);
    // for x < xlim do a the exp(-x)/exp(-x) multiply and a Taylor series expansion
    double lower_y = exp(-lower_x);
    double upper_y = exp(-upper_x);
    double add_lower , add_upper;
    if(lower_x >= 1.e-5)
     add_lower = NORM_FACTOR*lower_y*pow(lower_x,4)/(1.-lower_y);
    else 
     add_lower = NORM_FACTOR*pow(lower_x,3)*(1.- 0.5*lower_x); // T.S. + H.O.T.

    if(upper_x >= 1.e-5)
     add_upper = NORM_FACTOR*upper_y*pow(upper_x,4)/(1.-upper_y);
    else 
     add_upper = NORM_FACTOR*pow(upper_x,3)*(1.- 0.5*upper_x); // T.S.+ H.O.T.
    // one term taylor series for small x
    Ensure (integral >= 0.0 && integral <= 1.0);
    double PL = integral;
    double ROSL = PL - (add_upper-add_lower);

    return ROSL;
}
//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Rosseland spectrum over a frequency group.
 *
 * This function integrates the normalized Rosseland that is defined:
 * \f[
 *    r(x) = \frac{15}{4\pi^4} \frac{x^4 e^x}{(e^x - 1)^2}
 * \f]
 * where 
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 * and 
 * \f[
 *    R(\nu, T)d\nu = \frac{4 acT^3}{4\pi}r(x)dx
 * \f]
 * where \f$R(\nu, T)\f$ is the Rosseland and is defined
 * \f[
 *    R(\nu, T) = \frac{\partial B(\nu, T)}{\partial T}
 * \f]
 * \f[
 *    B(\nu, T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}} - 1}
 * \f]
 * The normalized Rosseland, integrated from 0 to \f$\infty\f$, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Rosseland function may integrate to something less than one.
 *
 * This function performs the following integration:
 * \f[
 *      \int_{x_1}^{x_N} r(x) dx
 * \f]
 * where \f$x_1\f$ is the low frequency bound and \f$x_N\f$ is the high
 * frequency bound of the multigroup data set.  This integration uses the
 * method of B. Clark (JCP (70)/2, 1987).  We use a 10-term Polylogarithmic
 * expansion for the normalized Planckian, except in the low-x limit, where
 * we use a 21-term Taylor series expansion.
 *
 * For the Rosseland we can relate the group interval integration to the
 * Planckian group interval integration, by equation 27 in B. Clark paper.
 * \f[
 *     \int_{x_1}^{x_N} r(x) dx = \int_{x_1}^{x_N} b(x) dx
 *     - \frac{15}{4\pi^4} \frac{x^4}{e^x - 1}
 * \f]
 * Therefore our Rosslenad group integration function can simply wrap the
 * Planckian function integratePlanckSpectrum(lowFreq, highFreq, T)
 * 
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 * \f[
 *      \int_{\nu_1}^{\nu_N} R(\nu,T) d\nu
 * \f]
 * then you must multiply by a factor of \f$\frac{4 acT^3}{4\pi}\f$ where a is
 * the radiation constant.  If you want to evaluate expressions like the
 * following:
 *\f[
 *      \int_{4\pi} \int_{\nu_1}^{\nu_N} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$ 4 acT^3 \f$.
 *
 * In the limit of \f$ T \rightarrow 0, r(T) \rightarrow 0\f$, therefore we
 * return a hard zero for a temperature equal to a hard zero.
 
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return integrated normalized Rosseland from x_low to x_high
 *
 *  
 * If groupindex =0  then an exception is thrown.
 *
 * \param groupIndex index of the frequency group to integrate [1,num_groups]
 * \param T          the temperature in keV (must be greater than 0.0)
 * \return           integrated normalized Plankian over the group specified
 *                   by groupIndex.
 *

 */

double CDI::integrateRosselandSpectrum(const int groupIndex,const double T)
{
    Insist  (!frequencyGroupBoundaries.empty(), "No groups defined!");
    Require (T >= 0.0);
    Require (groupIndex > 0 && 
	     groupIndex <= frequencyGroupBoundaries.size() - 1);

    // first determine the group boundaries for groupIndex
    double lowFreq = frequencyGroupBoundaries[groupIndex-1];
    double highFreq = frequencyGroupBoundaries[groupIndex];
    Check (highFreq  > lowFreq);

    // determine the upper and lower x
    double lower_x = lowFreq  / T;
    double upper_x = highFreq / T;
    double NORM_FACTOR = 0.25*coeff; //  15./(4.*PI4);

    // calculate the Planckian integral 
    double integral = integratePlanckSpectrum(lowFreq, highFreq, T);
    
    // Calculate the addition to the Planck integral 
    // for x > xlim=1.e-5
    //  double add_lower = NORM_FACTOR*pow(lower_x,4)/(exp(lower_x)-1.);
    //  double add_upper = NORM_FACTOR*pow(upper_x,4)/(exp(upper_x)-1.);
    // for x < xlim do a the exp(-x)/exp(-x) multiply and a Taylor series expansion
    double lower_y = exp(-lower_x);
    double upper_y = exp(-upper_x);
    double add_lower , add_upper;
    if(lower_x >= 1.e-5)
     add_lower = NORM_FACTOR*lower_y*pow(lower_x,4)/(1.-lower_y);
    else 
     add_lower = NORM_FACTOR*pow(lower_x,3)*(1.- 0.5*lower_x); // T.S. + H.O.T.

    if(upper_x >= 1.e-5)
     add_upper = NORM_FACTOR*upper_y*pow(upper_x,4)/(1.-upper_y);
    else 
     add_upper = NORM_FACTOR*pow(upper_x,3)*(1.- 0.5*upper_x); // T.S.+ H.O.T.
    // one term taylor series for small x
    Ensure (integral >= 0.0 && integral <= 1.0);
    double PL = integral;
    double ROSL = PL - (add_upper-add_lower);

    return ROSL;
}

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Planckian and Rosseland spectrum over a frequency group.
 *
 * This function integrates the normalized Rosseland that is defined:
 * \f[
 *    r(x) = \frac{15}{4\pi^4} \frac{x^4 e^x}{(e^x - 1)^2}
 * \f]
 * where 
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 * and 
 * \f[
 *    R(\nu, T)d\nu = \frac{4 acT^3}{4\pi}r(x)dx
 * \f]
 * where \f$R(\nu, T)\f$ is the Rosseland and is defined
 * \f[
 *    R(\nu, T) = \frac{\partial B(\nu, T)}{\partial T}
 * \f]
 * \f[
 *    B(\nu, T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}} - 1}
 * \f]
 * The normalized Rosseland, integrated from 0 to \f$\infty\f$, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Rosseland function may integrate to something less than one.
 *
 * This function performs the following integration:
 * \f[
 *      \int_{x_1}^{x_N} r(x) dx
 * \f]
 * where \f$x_1\f$ is the low frequency bound and \f$x_N\f$ is the high
 * frequency bound of the multigroup data set.  This integration uses the
 * method of B. Clark (JCP (70)/2, 1987).  We use a 10-term Polylogarithmic
 * expansion for the normalized Planckian, except in the low-x limit, where
 * we use a 21-term Taylor series expansion.
 *
 * For the Rosseland we can relate the group interval integration to the
 * Planckian group interval integration, by equation 27 in B. Clark paper.
 * \f[
 *     \int_{x_1}^{x_N} r(x) dx = \int_{x_1}^{x_N} b(x) dx
 *     - \frac{15}{4\pi^4} \frac{x^4}{e^x - 1}
 * \f]
 * Therefore our Rosslenad group integration function can simply wrap the
 * Planckian function integratePlanckSpectrum(lowFreq, highFreq, T)
 * 
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 * \f[
 *      \int_{\nu_1}^{\nu_N} R(\nu,T) d\nu
 * \f]
 * then you must multiply by a factor of \f$\frac{4 acT^3}{4\pi}\f$ where a is
 * the radiation constant.  If you want to evaluate expressions like the
 * following:
 *\f[
 *      \int_{4\pi} \int_{\nu_1}^{\nu_N} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$ 4 acT^3 \f$.
 *
 * In the limit of \f$ T \rightarrow 0, r(T) \rightarrow 0\f$, therefore we
 * return a hard zero for a temperature equal to a hard zero.
 
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return void the integrated normalized Planckian and Rosseland from x_low to x_high
 * are passed by reference in the function call
 *
 */

void CDI::integrateRosselandSpectrum(const double lowFreq,const double highFreq,const double T, double& PL, double& ROSL)
{
    Require (lowFreq >= 0.0);
    Require (highFreq >= lowFreq);
    Require (T >= 0.0);
    // return 0 if temperature is a hard zero
    if (T == 0.0)
    {
        PL = 0.; ROSL =0.;
	return ;
    }

    // determine the upper and lower x
    double lower_x = lowFreq  / T;
    double upper_x = highFreq / T;
    double NORM_FACTOR = 0.25*coeff; //  15./(4.*PI4);

    // calculate the Planckian integral 
    double integral = integratePlanckSpectrum(lowFreq, highFreq, T);
    
    // Calculate the addition to the Planck integral 
    // for x > xlim
    //  double add_lower = NORM_FACTOR*pow(lower_x,4)/(exp(lower_x)-1.);
    //  double add_upper = NORM_FACTOR*pow(upper_x,4)/(exp(upper_x)-1.);
    // for x < xlim do a exp(-x)/exp(-x) multiply
    double lower_y = exp(-lower_x);
    double upper_y = exp(-upper_x);
    double add_lower , add_upper;
    if(lower_x >= 1.e-5)
     add_lower = NORM_FACTOR*lower_y*pow(lower_x,4)/(1.-lower_y);
    else 
     add_lower = NORM_FACTOR*pow(lower_x,3)*(1.- 0.5*lower_x); // T.S. + H.O.T.

    if(upper_x >= 1.e-5)
     add_upper = NORM_FACTOR*upper_y*pow(upper_x,4)/(1.-upper_y);
    else 
     add_upper = NORM_FACTOR*pow(upper_x,3)*(1.- 0.5*upper_x); // T.S.+ H.O.T.
    // one term taylor series for small x
    Ensure (integral >= 0.0 && integral <= 1.0);
    PL = integral;
    ROSL = PL - (add_upper-add_lower);
    
    }

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Planckian and Rosseland spectrum over a frequency group.
 *
 * This function integrates the normalized Rosseland that is defined:
 * \f[
 *    r(x) = \frac{15}{4\pi^4} \frac{x^4 e^x}{(e^x - 1)^2}
 * \f]
 * where 
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 * and 
 * \f[
 *    R(\nu, T)d\nu = \frac{4 acT^3}{4\pi}r(x)dx
 * \f]
 * where \f$R(\nu, T)\f$ is the Rosseland and is defined
 * \f[
 *    R(\nu, T) = \frac{\partial B(\nu, T)}{\partial T}
 * \f]
 * \f[
 *    B(\nu, T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}} - 1}
 * \f]
 * The normalized Rosseland, integrated from 0 to \f$\infty\f$, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Rosseland function may integrate to something less than one.
 *
 * This function performs the following integration:
 * \f[
 *      \int_{x_1}^{x_N} r(x) dx
 * \f]
 * where \f$x_1\f$ is the low frequency bound and \f$x_N\f$ is the high
 * frequency bound of the multigroup data set.  This integration uses the
 * method of B. Clark (JCP (70)/2, 1987).  We use a 10-term Polylogarithmic
 * expansion for the normalized Planckian, except in the low-x limit, where
 * we use a 21-term Taylor series expansion.
 *
 * For the Rosseland we can relate the group interval integration to the
 * Planckian group interval integration, by equation 27 in B. Clark paper.
 * \f[
 *     \int_{x_1}^{x_N} r(x) dx = \int_{x_1}^{x_N} b(x) dx
 *     - \frac{15}{4\pi^4} \frac{x^4}{e^x - 1}
 * \f]
 * Therefore our Rosslenad group integration function can simply wrap the
 * Planckian function integratePlanckSpectrum(lowFreq, highFreq, T)
 * 
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 * \f[
 *      \int_{\nu_1}^{\nu_N} R(\nu,T) d\nu
 * \f]
 * then you must multiply by a factor of \f$\frac{4 acT^3}{4\pi}\f$ where a is
 * the radiation constant.  If you want to evaluate expressions like the
 * following:
 *\f[
 *      \int_{4\pi} \int_{\nu_1}^{\nu_N} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$ 4 acT^3 \f$.
 *
 * In the limit of \f$ T \rightarrow 0, r(T) \rightarrow 0\f$, therefore we
 * return a hard zero for a temperature equal to a hard zero.
 
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return void the integrated normalized Planckian and Rosseland from x_low to x_high
 * are passed by reference in the function call
 * * 
 * if  groupIndex=0  an exception is thrown.
 *
 * \param groupIndex index of the frequency group to integrate [1,num_groups]
 * \param T          the temperature in keV (must be greater than 0.0)
 * \return           integrated normalized Plankian over the group specified
 *                   by groupIndex.
 *

 */

void CDI::integrate_Rosseland_Planckian_Spectrum(const int groupIndex,const double T, double& PL, double& ROSL)
{
    Insist  (!frequencyGroupBoundaries.empty(), "No groups defined!");
    Require (T >= 0.0);
    Require (groupIndex > 0 && 
	     groupIndex <= frequencyGroupBoundaries.size() - 1);

    // first determine the group boundaries for groupIndex
    double lowFreq = frequencyGroupBoundaries[groupIndex-1];
    double highFreq = frequencyGroupBoundaries[groupIndex];
    Check (highFreq  > lowFreq);
 
    // return 0 if temperature is a hard zero
    if (T == 0.0)
    {
        PL = 0.; ROSL =0.;
	return ;
    }

    // determine the upper and lower x
    double lower_x = lowFreq  / T;
    double upper_x = highFreq / T;
    double NORM_FACTOR = 0.25*coeff; //  15./(4.*PI4);

    // calculate the Planckian integral 
    double integral = integratePlanckSpectrum(lowFreq, highFreq, T);
    
    // Calculate the addition to the Planck integral 
    // for x > xlim
    //  double add_lower = NORM_FACTOR*pow(lower_x,4)/(exp(lower_x)-1.);
    //  double add_upper = NORM_FACTOR*pow(upper_x,4)/(exp(upper_x)-1.);
    // for x < xlim do a exp(-x)/exp(-x) multiply
    double lower_y = exp(-lower_x);
    double upper_y = exp(-upper_x);
    double add_lower , add_upper;
    if(lower_x >= 1.e-5)
     add_lower = NORM_FACTOR*lower_y*pow(lower_x,4)/(1.-lower_y);
    else 
     add_lower = NORM_FACTOR*pow(lower_x,3)*(1.- 0.5*lower_x); // T.S. + H.O.T.

    if(upper_x >= 1.e-5)
     add_upper = NORM_FACTOR*upper_y*pow(upper_x,4)/(1.-upper_y);
    else 
     add_upper = NORM_FACTOR*pow(upper_x,3)*(1.- 0.5*upper_x); // T.S.+ H.O.T.
    // one term taylor series for small x
    Ensure (integral >= 0.0 && integral <= 1.0);
    PL = integral;
    ROSL = PL - (add_upper-add_lower);
    
    }

//---------------------------------------------------------------------------//
/*!
 *
 * \brief Integrate the Planckian and Rosseland spectrum over a frequency group.
 *
 * This function integrates the normalized Rosseland that is defined:
 * \f[
 *    r(x) = \frac{15}{4\pi^4} \frac{x^4 e^x}{(e^x - 1)^2}
 * \f]
 * where 
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 * and 
 * \f[
 *    R(\nu, T)d\nu = \frac{4 acT^3}{4\pi}r(x)dx
 * \f]
 * where \f$R(\nu, T)\f$ is the Rosseland and is defined
 * \f[
 *    R(\nu, T) = \frac{\partial B(\nu, T)}{\partial T}
 * \f]
 * \f[
 *    B(\nu, T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}} - 1}
 * \f]
 * The normalized Rosseland, integrated from 0 to \f$\infty\f$, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Rosseland function may integrate to something less than one.
 *
 * This function performs the following integration:
 * \f[
 *      \int_{x_1}^{x_N} r(x) dx
 * \f]
 * where \f$x_1\f$ is the low frequency bound and \f$x_N\f$ is the high
 * frequency bound of the multigroup data set.  This integration uses the
 * method of B. Clark (JCP (70)/2, 1987).  We use a 10-term Polylogarithmic
 * expansion for the normalized Planckian, except in the low-x limit, where
 * we use a 21-term Taylor series expansion.
 *
 * For the Rosseland we can relate the group interval integration to the
 * Planckian group interval integration, by equation 27 in B. Clark paper.
 * \f[
 *     \int_{x_1}^{x_N} r(x) dx = \int_{x_1}^{x_N} b(x) dx
 *     - \frac{15}{4\pi^4} \frac{x^4}{e^x - 1}
 * \f]
 * Therefore our Rosslenad group integration function can simply wrap the
 * Planckian function integratePlanckSpectrum(lowFreq, highFreq, T)
 * 
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 * \f[
 *      \int_{\nu_1}^{\nu_N} R(\nu,T) d\nu
 * \f]
 * then you must multiply by a factor of \f$\frac{4 acT^3}{4\pi}\f$ where a is
 * the radiation constant.  If you want to evaluate expressions like the
 * following:
 *\f[
 *      \int_{4\pi} \int_{\nu_1}^{\nu_N} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$ 4 acT^3 \f$.
 *
 * In the limit of \f$ T \rightarrow 0, r(T) \rightarrow 0\f$, therefore we
 * return a hard zero for a temperature equal to a hard zero.
 
 * \param T the temperature in keV (must be greater than 0.0)
 * 
 * \return void the integrated normalized Planckian and Rosseland from x_low to x_high
 * are passed by reference in the function call
 *
 */

void CDI::integrate_Rosseland_Planckian_Spectrum(const double lowFreq,const double highFreq,const double T, double& PL, double& ROSL)
{
    Require (lowFreq >= 0.0);
    Require (highFreq >= lowFreq);
    Require (T >= 0.0);
    // return 0 if temperature is a hard zero
    if (T == 0.0)
    {
        PL = 0.; ROSL =0.;
	return ;
    }

    // determine the upper and lower x
    double lower_x = lowFreq  / T;
    double upper_x = highFreq / T;
    double NORM_FACTOR = 0.25*coeff; //  15./(4.*PI4);

    // calculate the Planckian integral 
    double integral = integratePlanckSpectrum(lowFreq, highFreq, T);
    
    // Calculate the addition to the Planck integral 
    // for x > xlim
    //  double add_lower = NORM_FACTOR*pow(lower_x,4)/(exp(lower_x)-1.);
    //  double add_upper = NORM_FACTOR*pow(upper_x,4)/(exp(upper_x)-1.);
    // for x < xlim do a exp(-x)/exp(-x) multiply
    double lower_y = exp(-lower_x);
    double upper_y = exp(-upper_x);
    double add_lower , add_upper;
    if(lower_x >= 1.e-5)
     add_lower = NORM_FACTOR*lower_y*pow(lower_x,4)/(1.-lower_y);
    else 
     add_lower = NORM_FACTOR*pow(lower_x,3)*(1.- 0.5*lower_x); // T.S. + H.O.T.

    if(upper_x >= 1.e-5)
     add_upper = NORM_FACTOR*upper_y*pow(upper_x,4)/(1.-upper_y);
    else 
     add_upper = NORM_FACTOR*pow(upper_x,3)*(1.- 0.5*upper_x); // T.S.+ H.O.T.
    // one term taylor series for small x
    Ensure (integral >= 0.0 && integral <= 1.0);
    PL = integral;
    ROSL = PL - (add_upper-add_lower);
    
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
