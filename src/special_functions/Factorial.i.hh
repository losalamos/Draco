//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/factorial.i.hh
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide implementation of templatized factorial function.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef sf_factorial_i_hh
#define sf_factorial_i_hh

#include "ds++/Assert.hh"
#include "Factorial.hh"

namespace rtt_sf
{

//---------------------------------------------------------------------------//
/*! 
 * \brief factorial
 * 
 * This is a recursive function that computes k!
 *
 * \param k A number whose factorial is desired.
 * \return \f$n!\f$
 * \post \c Result>=1
 */
template< typename T >
T factorial( T const k ) 
{
    // only initialize this once (keyword: static)
    static T const tabularValue[] =
    {
        1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880,
        3628800, 39916800, 479001600
        // These are the next two values in the series.  Unfortunately, they
        // are too big to be held by type long.
        // , 6227020800, 87178291200
    };
    static unsigned const N(13);

    if( k <= 1 ) 
	return 1;
    if( k < N )
        return tabularValue[k];
    else
    {
        Check(k>=N);

    // recursive algorithm:
    // return  k * factorial(k-1);

    // direct algorithm
    
        T Result( tabularValue[N-1] );
        for( unsigned i=N; i<=k; i++)
            Result *= i;
        Ensure(Result>tabularValue[N-1]);
        return Result;
    }
}

} // end namespace rtt_sf

#endif // sf_factorial_i_hh

//---------------------------------------------------------------------------//
//              end of sf/factorial.i.hh
//---------------------------------------------------------------------------//
