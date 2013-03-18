//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/pythag.hh
 * \author Kent Budge
 * \date   Mon Aug  9 13:17:31 2004
 * \brief  Calculate hypotenuse of a right triangle.
 * \note   Copyright (C) 2004-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_pythag_hh
#define dsxx_pythag_hh

#include <cmath>
#include "square.hh"
#include "abs.hh"

namespace rtt_dsxx
{
//---------------------------------------------------------------------------//
/*! 
 * \brief Compute the hypotenuse of a right triangle.
 *
 * This function evaluates the expression \f$\sqrt{a^2+b^2}\f$ in a way that
 * is insensitive to roundoff and preserves range.
 *
 * \arg \a Real A real number type
 * \param a First leg of triangle
 * \param b Second leg of triangle
 * \return Hypotenuse, \f$\sqrt{a^2+b^2}\f$
 */

template<class Real>
inline double pythag(Real a, Real b)
{
    Real absa = abs(a), absb = abs(b);
    if (absa>absb)
    {
	return absa*std::sqrt(1.0+square(absb/absa));
    }
    else
    {
	if (absb==0.0)
	{
	    return 0.0;
	}
	else
	{
	    return absb*std::sqrt(1.0+square(absa/absb));
	}
    }
}

} // end namespace rtt_dsxx

#endif // dsxx_pythag_hh

//---------------------------------------------------------------------------//
// end of ds++/pythag.hh
//---------------------------------------------------------------------------//
