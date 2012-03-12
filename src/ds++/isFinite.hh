//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   ds++/isFinite.hh
 * \brief  Checks on floating point numbers
 *
 * This header defines several portable functions that check on the status of
 * floating point numbers.
 */
//---------------------------------------------------------------------------//
// $Id$ 
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_isFinite_hh
#define rtt_dsxx_isFinite_hh

#include "Assert.hh"
#include <algorithm>
#include <iterator>
#include <functional>

namespace rtt_dsxx
{
/*---------------------------------------------------------------------------*/
/*! 
 * \brief Return true if the value x is NaN.
 * 
 * \param x value to check for NaN-ness.
 * \return \c true if x is NaN; \c false otherwise.
 */
template< typename T >
inline bool isNaN( T const & x )
{
    return ( (x == x) == false );
}

//---------------------------------------------------------------------------//
//! Return true if x equals positive or negative infinity.
template< typename T > DLL_PUBLIC
bool isInfinity( T const & x );


/*---------------------------------------------------------------------------*/
/*!
 * \brief Return true if the value x is a finite number.
 * 
 * \param x value to check for finite-ness.
 * \return bool - \c true if x is finite.
 *
 * Taken from http://beagle.gel.ulaval.ca/refmanual/a01120.html
 */
template< typename T >
inline bool isFinite( T const & x )
{
    return ((isNaN<T>(x) == false ) && (isInfinity(x) == false ));
}

} // ane of namespace rtt_dsxx

#endif // rtt_dsxx_isFinite_hh
//---------------------------------------------------------------------------//
//                           end of isFinite.hh
//---------------------------------------------------------------------------//



