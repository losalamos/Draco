//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadServices.i.hh
 * \author Kelly Thompson
 * \date   Mon Nov  8 14:23:03 2004
 * \brief  Member definitions of class QuadServices
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_QuadServices_i_hh
#define quadrature_QuadServices_i_hh

#include "QuadServices.hh"

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*! 
 * \brief kronecker_delta
 * 
 * Return 1 if test_value == offset, otherwise return 0;
 * 
 * \param test_value
 * \param offset
 * \return 1 if test_value == offset, otherwise return 0;
 *
 * \bug move to ds++
 */
template< typename T >
T QuadServices::kronecker_delta( T const test_value,
				 T const offset ) const
{ 
    return (test_value==offset) ? 1 : 0; 
}

//---------------------------------------------------------------------------//
/*! 
 * \brief factorial
 * 
 * This is a recursive function that computes k!
 *
 * \param k 
 * \return k!
 *
 * \bug move to ds++
 */
template< typename T >
T QuadServices::factorial( T const k ) const
{
    if( k <= 1 ) 
	return 1;
    else 
	return  k * factorial(k-1);
}

} // end namespace rtt_quadrature

#endif // quadrature_QuadServices_i_hh

//---------------------------------------------------------------------------//
//              end of quadrature/QuadServices.i.hh
//---------------------------------------------------------------------------//
