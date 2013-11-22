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
#include "ds++/config.h" // defines DRACO_DIAGNOSTICS

namespace rtt_memory
{

bool set_memory_checking(bool new_status);

unsigned total_allocation();
unsigned peak_allocation();
unsigned largest_allocation();

} // namespace rtt_memory

#endif

//---------------------------------------------------------------------------//
// end of memory/memory.hh
//---------------------------------------------------------------------------//
