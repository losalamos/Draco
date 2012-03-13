//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   ds++/isFinite.i.hh
 * \brief  Template implementation for isFinite
 * \note   Copyright (C) 2006-2012 Los Alamos National Security, LLC 
 * \version $Id$
 */
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_isFinite_i_hh
#define rtt_dsxx_isFinite_i_hh

#include "Assert.hh"
#include "Soft_Equivalence.hh"
#include <algorithm>
#include <iterator>
#include <functional>
#include <limits>

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



