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

#include "rng/Names.hh"
#include "rng/Rnd_Control.hh"
#include "rng/Sprng.hh"
#include "sprng.h"

RNGSPACE

//---------------------------------------------------------------------------//
// CONSTANTS USED BY SPRNG RANDOM CLASS CLIENTS
//---------------------------------------------------------------------------//

const int max_buffer = MAX_PACKED_LENGTH;

//---------------------------------------------------------------------------//
// SPRNG functions set in this namespace
//---------------------------------------------------------------------------//

using ::pack_sprng;
using ::unpack_sprng;
using ::init_sprng;
using ::sprng;
using ::print_sprng;
using ::spawn_sprng;
using ::free_sprng;

CSPACE

#endif                          // __rng_Random_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Random.hh
//---------------------------------------------------------------------------//
