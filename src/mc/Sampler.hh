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

} // end namespace sampler
} // end namespace rtt_mc

#endif                          // __mc_Sampler_hh__

//---------------------------------------------------------------------------//
//                              end of mc/Sampler.hh
//---------------------------------------------------------------------------//
