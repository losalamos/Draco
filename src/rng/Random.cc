//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Random.cc
 * \author Thomas M. Evans
 * \date   Mon Jun 15 08:18:20 1998
 * \brief  Variable definitions in the rtt_rng namespace.
 */ 
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Random.hh"

namespace rtt_rng 
{

//---------------------------------------------------------------------------//
// VARIABLE DEFINITIONS IN RNG
//---------------------------------------------------------------------------//

/*!
 * \brief Random number stream index definition.
 *
 * This variable is a global-style index that can be used to reference the
 * current number of random number stream states.  For example, if 100 random
 * number stream states have been used across N processors, setting
 * rn_stream=100 on each processor informs all processes that 100 random
 * number states have been assigned.  Thus, reproducibility can be maintained 
 * across the entire problem domain.
 *
 * Remember, one can use any variable to perform this task.  We provide this
 * namespace-scoped variable as a convenience.
 */ 
int rn_stream = 0;

} // end namespace rtt_rng

//---------------------------------------------------------------------------//
//                              end of Random.cc
//---------------------------------------------------------------------------//
