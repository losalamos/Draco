//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Random.hh
 * \author Thomas M. Evans
 * \date   Tue May 12 09:00:52 1998
 * \brief  General header file for including the Rnd_Control, Sprng, and
 *         SPRNG library headers.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __rng_Random_hh__
#define __rng_Random_hh__

//===========================================================================//
//
// revision history:
// -----------------
// 0)  original
// 1)  1-SEP-99 : added doxygen comments
// 
//===========================================================================//

// get configurations for the RNG package
#include <rng/config.h>

// Loads the primary rng package files
#include "Rnd_Control.hh"
#include "Sprng.hh"

// Header file that includes the SPRNG library.
#include "rng_sprng.h"

namespace rtt_rng 
{

//---------------------------------------------------------------------------//
// CONSTANTS USED BY SPRNG RANDOM CLASS CLIENTS
//---------------------------------------------------------------------------//

//! Max size of a Sprng random number state.
const int max_buffer = MAX_PACKED_LENGTH;

//---------------------------------------------------------------------------//
// variables useful for RNG users
//---------------------------------------------------------------------------//

//! Random number stream index declaration.
extern int rn_stream;

} // end namespace rtt_rng

#endif                          // __rng_Random_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Random.hh
//---------------------------------------------------------------------------//
