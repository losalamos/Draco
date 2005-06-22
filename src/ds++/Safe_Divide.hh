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
     * This function requires an extra multiply to compute a limit on the
     * numerator, Note that limit*denominator may overflow on some
     * archetectures for denominatior > 1. 
     */

    FT limit = std::numeric_limits<FT>::max();
    FT result = (numerator < limit * denominator) ? numerator / denominator : limit;

    return result;
}

} // end namespace rtt_dsxx

#endif // dsxx_Save_Divide_hh

//---------------------------------------------------------------------------//
//              end of ds++/Save_Divide.hh
//---------------------------------------------------------------------------//
