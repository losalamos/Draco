//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Sampler.hh
 * \author Todd J. Urbatsch
 * \date   April 4, 2000 
 * \brief  Collection of sampling functions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __mc_Sampler_hh__
#define __mc_Sampler_hh__

//===========================================================================//
// namespace Sampler - 
//
// Purpose : holds sampling functions
//
// revision history:
// -----------------
// 0) 04-04-00 : original
// 1) 08-11-00 : modified sample_general_linear fn to allow for, in addition
//               to a totally positive function, a totally negative or a
//               zero-valued function.
//
//===========================================================================//

#include "ds++/Assert.hh"

namespace rtt_mc 
{
namespace sampler
{

//---------------------------------------------------------------------------//
// SAMPLE_GENERAL_LINEAR
//---------------------------------------------------------------------------//
/*!
 * \brief Sample x in [a,b] from a general linear pdf, f(x).
 *
 * The linear function is broken into two linear functions f(x) = n(x) +
 * p(x), one with negative slope and one with positive slope, and each with
 * zero intercept at b and a respectively.  f(x) does not have to be
 * normalized.  Uses two random numbers on each call.
 *
 * \param ran random number object
 * \param a lower extent of independent variable range.
 * \param b upper extent of independent variable range.
 * \param f_of_a f(a).
 * \param f_of_b f(b). 
 * \return x \f$ \in \f$ [a,b] sampled from \f$ f(x) \f$. 
 */
template<class Ran>
inline double sample_general_linear(Ran ran, const double a, const double b, 
				    const double f_of_a, const double f_of_b) 
{
    using namespace std;

    // require that a<b and the function is nonnegative (could require that
    // it is either nonnegative OR nonpositive, but we will stick with >=0).
    Require (a < b);
    Require ((f_of_a >= 0.0 && f_of_b >= 0.0) || 
	     (f_of_a <= 0.0 && f_of_b <= 0.0));

    // return value of x in [a,b], sampled from f(x).
    double x;

    // calculate the probability of selecting the function with negative
    // slope.
    double prob_neg;
    if (std::fabs(f_of_a + f_of_b) > 0)
	prob_neg = f_of_a / (f_of_a + f_of_b);
    else 
	prob_neg = 0.5;

    // calculate length of independent extent.
    double delta_x = b - a;

    // sample which function to further sample.
    if (ran.ran() <= prob_neg)
    {
	// sample x from function with negative slope, n(x).
	x = b - delta_x * sqrt(ran.ran());
    }
    else
    {
	// sample x from function with positive slope, p(x).
	x = a + delta_x * sqrt(ran.ran());
    }

    Ensure (x >= a && x <= b);

    return x;
}

//---------------------------------------------------------------------------//
// SAMPLE_PLANCKIAN_FREQUENCY
//---------------------------------------------------------------------------//
/*!
 * \brief Sample a frequency from the continuous Planckian distribution.
 *
 * Given a temperature kT, a frequency in the same units as kT is sampled
 * from a continuous Planckian distribution.  
 *
 * This function utilizes Barnett and Canfield's truncated infinite series
 * technique [UCIR-473 (UCRL-125393)].  They begin by writing the normalized
 * Planckian, b(x) = (15/pi^4) x^3/(e^x-1), as Sum_{n=1}^{\infty} p_n f_n(x),
 * where x=hnu/kT is the reduced frequency, p_n = 90/(pi^4 n^4), f_n(x) =
 * (n^4/6) x^3 e^{-nx} and Sum_{n=1}^infty p_n = 1.  Thus, with probability
 * p_n we can sample f_n(x) for x.  Sampling n from p_n requires satisfying
 * (pi^4/90)ran <= Sum_{k=1}^{n} 1/k^4.  The search usually only requires one
 * iteration (i.e., n=1), but, without machine-error, the search can be
 * unbounded.  The search is bounded by truncating (rounding down) the
 * decimal representation of (pi^4/90) = 1.0823232337....  The maximum number
 * of iterations goes as follows:
 *
 *     Significant Digits   Truncated $\pi^4/90$   Max iters
 *
 *                  6         1.08232              47   
 *                  7         1.082323             113  
 *                  8         1.0823232            215  
 *                  9         1.08232323           448  
 *                 10         1.082323233          775  
 *                 11         1.0823232337         2744 
 *
 * Large n usually corresponds to a small frequency, so the concern over
 * significant digits is largely moot.  However, given that the maximum
 * number of iterations is rarely encountered, the computational cost is not
 * excessive (but the expenditure of random numbers might get wasteful).  
 *
 * Once n is found, x is sampled from f_n(x) according to the procedure in
 * Everett and Cashwell's "A Third Monte Carlo Sampler":
 * 
 * x = -1/n * ln(ran1 * ran2 * ran3 * ran4)
 *
 * The frequency, then, is h*nu = x*kT.
 *
 * The random number object must have a member function called "ran()" that
 * returns a random number of type double between 0 and 1.
 *
 * \param ran random number object
 * \param k_temperature temperature (kT)
 * \return hnu frequency in same units as k_temperature 
 */
template<class Ran>
inline double sample_planckian_frequency(Ran ran, const double k_temperature)
{
    using namespace std;
    using rtt_mc::global::soft_equiv;

    // check that temp is nonnegative
    Require (k_temperature >= 0.0);

    // declaration of return value of frequency
    double hnu;

    // initialize iteration, conditional probability and its sum
    double n       = 1.0;
    double p_n     = 1.0;
    double sum_p_n = 1.0;
 
    // 8 significant digits in the search for n
    double sampled_cdf = 1.0823232 * ran.ran();
    Check (sampled_cdf > 0.0);
 
    // continue the truncated infinite series sum if necessary
    while (sampled_cdf > sum_p_n)
    {
        n       = n + 1.0;
        p_n     = 1/n;
        sum_p_n = sum_p_n + p_n*p_n*p_n*p_n;
    }

    // given n, sample nu (note that these are four different random numbers)
    hnu = -p_n * k_temperature * log(ran.ran()*ran.ran()*ran.ran()*ran.ran());  
 
    // nu should only be zero if kT=0 or if all four of the random numbers
    // above are 1.0; otherwise nu should be > 0.  nu>=0.0 is a loose check. 
    Ensure (hnu >= 0.0);

    // return the frequency in the same units as the temperature
    return hnu;
}

} // end namespace sampler
} // end namespace rtt_mc

#endif                          // __mc_Sampler_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Sampler.hh
//---------------------------------------------------------------------------//
