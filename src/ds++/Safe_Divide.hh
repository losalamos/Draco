//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Save_Divide.hh
 * \author Mike Buksas
 * \date   Tue Jun 21 15:35:05 2005
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
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

template <typename FT>
inline FT safe_pos_divide (const FT& numerator, const FT& denominator)
{

    /* Implement division which maxes out at std::numerics_limits<FT>::max()
     * when the demoninator is too small.
     *
     * The arguments are assumed to be positive
     *
     * Compared to straight division, this function requires an actra function
     * call to get the max value, and  extra multiply and comparison (via
     * std::max) to compute a limit on the numerator. The comparison prevents
     * numeric overflow limit * denominator, which would generate a
     * floating-point exception and really hose the performance.
     *
     * The result is only slightly more operations than a straight divide and
     * more readable than including if-blocks around all problematic divisions.
     *
     */

    FT limit = std::numeric_limits<FT>::max();
    FT numerator_bound = limit * std::min (1.0, denominator);
    FT result = (numerator < numerator_bound) ? numerator / denominator : limit;

    return result;
}

} // end namespace rtt_dsxx

#endif // dsxx_Save_Divide_hh

//---------------------------------------------------------------------------//
//              end of ds++/Save_Divide.hh
//---------------------------------------------------------------------------//
