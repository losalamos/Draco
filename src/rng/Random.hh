//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Random.hh
 * \author Thomas M. Evans
 * \date   Tue May 12 09:00:52 1998
 * \brief  General header file for including the Rnd_Control, Sprng, and
 *         SPRNG library headers.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_rng_Random_hh
#define rtt_rng_Random_hh

// get configurations for the RNG package
#include <rng/config.h>

// Header file that includes the SPRNG library.
#include "rng_sprng.h"

// Loads the primary rng package files
#include "Rnd_Control.hh"
#include "Sprng.hh"

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

#endif                          // rtt_rng_Random_hh

//---------------------------------------------------------------------------//
//                              end of rng/Random.hh
//---------------------------------------------------------------------------//
