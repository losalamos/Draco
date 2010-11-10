//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/sign.hh
 * \author Kent Budge
 * \date   Wed Jul  7 09:14:09 2004
 * \brief  Reproduce the Fortran SIGN function.
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_sign_hh
#define dsxx_sign_hh

#include "abs.hh"

namespace rtt_dsxx
{
//-------------------------------------------------------------------------//
/*!
 * \brief  Transfer the sign of the second argument to the first argument.
 *
 * This is a replacement for the FORTRAN SIGN function.  It is
 * useful in numerical algorithms that are roundoff or overflow
 * insensitive and should not be deprecated.
 *
 * \arg \a Ordered_Group
 * A type for which \c operator< and unary \c operator- are defined and which
 * can be compared to literal \c 0.
 *
 * \param a
 * Argument supplying magnitude of result.
 *
 * \param b
 * Argument supplying sign of result.
 *
 * \return \f$|a|sgn(b)\f$
 */
 
template <class Ordered_Group>
inline Ordered_Group sign(Ordered_Group a,
			  Ordered_Group b)
{
    using rtt_dsxx::abs; // just to be clear
    
    if (b<0) 
    {
	return -abs(a);
    }
    else 
    {
	return abs(a);
    }
}


} // end namespace rtt_dsxx

#endif // dsxx_sign_hh

//---------------------------------------------------------------------------//
//              end of ds++/sign.hh
//---------------------------------------------------------------------------//
