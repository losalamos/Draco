//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   ds++/abs.hh
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Absolute value template
 *
 * Absolute values are a mess in the STL, in part because they are a mess
 * in the standard C library. We do our best to give a templatized version here.
 */
//---------------------------------------------------------------------------//
// $Id$ 
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_abs_hh
#define rtt_dsxx_abs_hh

#include <cmath>
#include <cstdlib>

namespace rtt_dsxx
{

//-------------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date   Thu Jan 23 08:41:54 MST 2003
 * \brief  Return the absolute value of the argument.
 * 
 * \arg \a Ordered_Group
 * A type for which operator< and unary operator- are defined.
 *
 * \param a
 * Argument whose absolute value is to be calculated.
 *
 * \return \f$|a|\f$
 */

template <class Ordered_Group>
inline Ordered_Group abs(Ordered_Group a)
{ 
    if (a<0)
    {
	return -a;
    }
    else
    {
	return a;
    }
}

//---------------------------------------------------------------------------//
/*!
 * Specializations for standard types
 */

/* There is no standard abs function for float -- one reason why we define
 * this template! */

/* double */

template <>
inline double abs(double a)
{
    return std::fabs(a);
}

/* For int */

template<>
inline int abs(int a)
{
    return std::abs(a);
}

/* For long */

template<>
inline long abs(long a)
{
    return std::labs(a);
}

}  // rtt_dsxx

#endif // rtt_dsxx_abs_hh
//---------------------------------------------------------------------------//
//                           end of abs.hh
//---------------------------------------------------------------------------//



