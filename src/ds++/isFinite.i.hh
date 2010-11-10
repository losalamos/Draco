//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   ds++/isFinite.i.hh
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Template implementation for isFinite
 */
//---------------------------------------------------------------------------//
// $Id$ 
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_isFinite_i_hh
#define rtt_dsxx_isFinite_i_hh

#include <algorithm>
#include <iterator>
#include <functional>
#include <limits>
#include "Assert.hh"
#include "Soft_Equivalence.hh"

namespace rtt_dsxx
{
/*---------------------------------------------------------------------------*/
/*! 
 * \brief Return true if x equals positive or negative infinity.
 * 
 * \param x value to check for infinite value.
 * \return \c true if x equals positive or negative infinity; \c false
 * otherwise.
 */
template< typename T >
bool isInfinity( T const & x )
{
    if( std::numeric_limits<T>::has_infinity)
    {
	return ( ( x == std::numeric_limits<T>::infinity() ) ||
		 ( x == -std::numeric_limits<T>::infinity() ) );
    }
    else
    {
        T lZero(0.0);
        T lInfinity(1.0/lZero);
        return (( x == lInfinity ) || ( x == -lInfinity ));
    }
}

}  // rtt_dsxx

#endif // rtt_dsxx_isFinite_i_hh
//---------------------------------------------------------------------------//
//                           end of isFinite.i.hh
//---------------------------------------------------------------------------//



