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
// 2) 12-19-01 : added function sample_bin_from_discrete_cdf and its
//               supporting function is_this_cdf_valid for DBC use.
//
//===========================================================================//

#include <vector>
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

namespace rtt_mc 
{
namespace sampler
{

// STL typedefs
typedef std::vector<double> sf_double;


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
 *
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

//---------------------------------------------------------------------------//
// IS_THIS_CDF_VALID
//---------------------------------------------------------------------------//
/*!
 * \brief Check if a discrete cdf (vector<double>) is valid.
 *
 * A discrete cumulative distribution function (cdf) is valid if each element
 * is nonnegative and not less than the previous value.  These criteria stem
 * from the criterion of probability density functions (pdf) that p(i)>0.  No
 * check is made on whether the cdf is normalized or not, because, although
 * it is technically required, it is oftentimes more numerically efficient to
 * sample an unnormalized cdf.
 *
 * \param cdf vector<double> cumulative distribution function, assumed to be
 * in increasing order
 *
 * \return bool whether discrete cdf is valid
 *
 * Limited to vector<double> now.  Add capability for other containers as
 * needed.
 *
 */
inline bool is_this_cdf_valid(const sf_double &cdf) 
{
    // make sure there is actually some data
    if (cdf.size() <= 0) return false;

    // check for negativities in the first value.
    if (cdf[0] < 0.0) return false;

    // continuing with the remaining values, check for negativities and
    // negative slopes in the cdf
    for (int i = 1; i < cdf.size(); i++)
	if (cdf[i] < 0.0 || cdf[i] < cdf[i-1]) return false;

    // if none of the checks fail, return true
    return true;
}

//---------------------------------------------------------------------------//
// SAMPLE_BIN_FROM_DISCRETE_CDF
//---------------------------------------------------------------------------//
/*!  
 * \brief Sample a bin, in [0, num_bins-1], from a discrete histogram
 *        cumulative distribution function (cdf). 
 *
 * \param ran random number object that has ran() member function
 * \param cdf vector<double> cumulative distribution function, constant in
 *            each bin  
 *
 * \return sampled bin, or index, of cdf in [0, num_bins-1]
 *
 * Limited to vector<double> now.  Add capability for other containers as
 * needed.
 *
 */
template<class Ran>
inline int sample_bin_from_discrete_cdf(Ran ran_gen, const sf_double &cdf) 
{
    using rtt_dsxx::soft_equiv;

    // check that the cdf is a valid cdf
    Check (is_this_cdf_valid(cdf));

    // locally set the number of bins
    int num_bins = cdf.size();

    // get a random number and multiply it by the cdf's normalization factor
    // if the cdf is not normalized.
    double xi = ran_gen.ran();

    if (!soft_equiv(cdf[num_bins-1], 1.0, 1.0e-12))
	xi *= cdf[num_bins-1];

    // perform a binary search for the bin, where 
    //    bin is the running index,
    //    bin_lo is the lower, inclusive bound, and 
    //    bin_hi is the upper, inclusive bound.
    // the initial value of bin is the middle bin for odd num_bins or the bin
    // just below the midpoint for even num_bins.
    int bin_lo = 0;
    int bin_hi = num_bins - 1;
    int bin    = bin_hi * 0.5;
    
    while (bin_lo != bin_hi)
    {
	Check (bin >= 0 && bin < num_bins);

	// if the random variable is below this cdf value, move the upper
	// (inclusive) bound to here.
	if (xi <= cdf[bin])
	    bin_hi = bin;

	// if the random variable is above this cdf value, move the lower
	// (inclusive) bound to the next bin up.
	else
	    bin_lo = bin + 1;

	// split the difference to get the next estimate of the bin
	bin = (bin_lo + bin_hi) * 0.5;
    }

    // make rudimentary checks on the bin and return
    Ensure (bin == bin_lo && bin == bin_hi);
    Ensure (bin >= 0      && bin < num_bins);

    return bin;
}

} // end namespace sampler
} // end namespace rtt_mc

#endif                          // __mc_Sampler_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Sampler.hh
//---------------------------------------------------------------------------//
