//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Save_Divide.hh
 * \author Mike Buksas
 * \date   Tue Jun 21 15:35:05 2005
 * \brief  Provide protected division functions.
 * \note   Copyright 2004 The Regents of the University of California.
 *
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_Save_Divide_hh
#define dsxx_Save_Divide_hh

#include "Soft_Equivalence.hh"
#include <limits>

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/**
 * \brief Provide a non-overflowing division routine
 * \arg dividend The number *being divided*
 * \arg divisor  the number *by which* the dividend is being divided
 * \return The quotient.
 *
 * Implement division  which maxes out at std::numerics_limits<FT>::max() when
 * the divisor is too small.
 *
 * The arguments are assumed to be positive
 *
 * The argument names are as follows: quotient = dividend / divisor. The
 * function is written so that a/b can be replaced with safe_pos_divide(a,b)
 * with the arguments in the same order.
 *
 * Compared to straight division, this function requires an extra multiply and
 * two comparisona The limit value for the data type must also be determined,
 * but this is only done once for each data type. 
 *
 * The result is only slightly more operations than a straight divide and more
 * readable than including if-blocks around all problematic divisions.
 * 
 */
template <typename FT>
inline FT safe_pos_divide (const FT& dividend, const FT& divisor)
{

    static FT limit = std::numeric_limits<FT>::max();
    const  FT dividend_bound = limit * std::min (1.0, divisor);
    return (dividend < dividend_bound) ? dividend / divisor : limit;

}

} // end namespace rtt_dsxx

#endif // dsxx_Save_Divide_hh

//---------------------------------------------------------------------------//
//              end of ds++/Save_Divide.hh
//---------------------------------------------------------------------------//
