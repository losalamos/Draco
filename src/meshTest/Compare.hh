//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/Compare.hh
 * \author Randy M. Roberts
 * \date   Wed Nov 24 09:40:43 1999
 * \brief  Compare inexact numeric results.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_Compare_hh__
#define __meshTest_Compare_hh__

#include <functional>
#include <limits>
#include <cmath>

namespace rtt_meshTest
{

//===========================================================================//
/*!
 * \class NearlyEqualToIsExact
 *
 * \brief Base class for inexact/exact comparison.
 *
 * The NearlyEqualToIsExact class is templated the type and a bool.
 * The bool is used to determine whether the comparison should be
 * exact or inexact.  This is accomplished through partial specialization.
 * This class should not be used directly, but through its derived class,
 * NearlyEqualTo<T>
 *
 * The constructor can take the number of digits of desired accuracy.
 * For the inexact specialization, this defaults to some value,
 * that at one time was 6, but may have changed.
 */
//===========================================================================//

template<bool IS_EXACT, class T>
struct NearlyEqualToIsExact
    : public std::binary_function<T,T,bool>
{
    NearlyEqualToIsExact() { /* empty */ }
    NearlyEqualToIsExact(int) { /* empty */ }

    bool operator()(T lhs, T rhs) const
    {
	return lhs == rhs;
    }
};

template<class T>
struct NearlyEqualToIsExact<false, T>
    : public std::binary_function<T,T,bool>
{
    int ndigits;

    static const int DefaultAccuracy = 6;
    
    NearlyEqualToIsExact(int ndigits_in = DefaultAccuracy)
	: ndigits(ndigits_in)
    {
	// empty
    }

    bool operator()(T lhs, T rhs) const
    {
	T zz;
    
	if (rhs != 0.)
	    zz = std::log10(std::abs((lhs - rhs)/rhs));
	else
	    zz = std::log10(std::abs(lhs));

	return zz <= -T(ndigits);
    }
};

//===========================================================================//
/*!
 * \class NearlyEqualTo
 *
 * \brief Class for inexact/exact comparison.
 *
 * The NearlyEqualTo<T> class is templated on the type of the comparison,
 * and is inherited from NearlyEqualToIsExact<bool,T>.
 * The bool is set to std::numeric_limits<T>::is_exact in order
 * to determine whether the comparison should be
 * exact or inexact.  This is accomplished through partial specialization
 * of the base class.
 *
 * The constructor can take the number of digits of desired accuracy.
 */
//===========================================================================//

template<class T>
struct NearlyEqualTo
    : public NearlyEqualToIsExact<std::numeric_limits<T>::is_exact,T>
{
    NearlyEqualTo()
    {
	// empty
    }
    NearlyEqualTo(int ndigits_in)
	: NearlyEqualToIsExact<std::numeric_limits<T>::is_exact,T>(ndigits_in)
    {
	// empty
    }
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_Compare_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/Compare.hh
//---------------------------------------------------------------------------//
