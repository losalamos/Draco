//----------------------------------*-C++-*----------------------------------//
// Random.hh
// Thomas M. Evans
// Tue May 12 09:00:52 1998
//---------------------------------------------------------------------------//
// @> General header file for including the SPRNG library classes
//---------------------------------------------------------------------------//

#ifndef __rng_Random_hh__
#define __rng_Random_hh__

//===========================================================================//
// Random - 
//
// Purpose : header file for including the SPRNG library classes in IMCTEST
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

// get configurations for the RNG package
#include <rng/config.h>

// get rng package files
#include "Rnd_Control.hh"
#include "Sprng.hh"

// Header file that includes the SPRNG library.
#include "rng_sprng.h"

namespace rtt_rng 
{

//---------------------------------------------------------------------------//
// CONSTANTS USED BY SPRNG RANDOM CLASS CLIENTS
//---------------------------------------------------------------------------//

const int max_buffer = MAX_PACKED_LENGTH;

//---------------------------------------------------------------------------//
// variables useful for RNG users
//---------------------------------------------------------------------------//

extern int rn_stream;

} // end namespace rtt_rng

#endif                          // __rng_Random_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Random.hh
//---------------------------------------------------------------------------//
