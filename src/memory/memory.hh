//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   memory/memory.hh
 * \author Kent G. Budge
 * \brief  Memory utilities for diagnostic purposes.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved.
 *
 * The memory utilities were written to address a need to identify the memory
 * "high-water mark" in a call sequence. This was not available with the
 * existing memory checking tools. Other capabilities gradually accreted
 * themselves to this set of utilities, such as leak characterization.
 */
//----------------------------------------------------------------------------//

#ifndef memory_memory_hh
#define memory_memory_hh

#include "ds++/config.h" // defines DRACO_DIAGNOSTICS
#include <cstdint>       // cstdint not available on PGI
#include <cstdlib>
#include <iostream>

namespace rtt_memory {

bool set_memory_checking(bool new_status);
uint64_t total_allocation();
uint64_t peak_allocation();
uint64_t largest_allocation();
void report_leaks(std::ostream &);

//! Register rtt_dsxx::print_stacktrace() as the respose to std::bad_alloc.
void out_of_memory_handler(void);

} // namespace rtt_memory

#endif

//----------------------------------------------------------------------------//
// end of memory/memory.hh
//----------------------------------------------------------------------------//
