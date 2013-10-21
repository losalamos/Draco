//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   memory/memory.hh
 * \author Kent G. Budge
 * \brief  memory utilities for diagnostic purposes
 * \note   Copyright (C) 2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: memory.hh 7239 2013-10-07 20:29:39Z kellyt $
//---------------------------------------------------------------------------//

#ifndef memory_memory_hh
#define memory_memory_hh

#include <cstdlib>

#include "diagnostics/config.h"

namespace rtt_memory
{

#if DRACO_DIAGNOSTICS & 2

bool set_memory_checking(bool new_status);

unsigned total_allocation();
unsigned peak_allocation();

#else

// If not wanted, the diagnostics are defined as empty inline functions so
// that the library, with the intercepting implementation of operator new,
// need not be linked.

inline bool set_memory_checking(bool){ return false; }

inline unsigned total_allocation(){ return 0; }
inline unsigned peak_allocation(){ return 0; }

#endif

} // namespace rtt_memory

#endif

//---------------------------------------------------------------------------//
// end of memory/memory.hh
//---------------------------------------------------------------------------//
