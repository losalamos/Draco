//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   memory/memory.cc
 * \author Kent G. Budge
 * \brief  memory diagnostic utilities
 * \note   Copyright (C) 2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: memory.cc 7133 2013-06-11 17:54:11Z kellyt $
//---------------------------------------------------------------------------//

#include <map>

#include "memory.hh"

#ifndef _GLIBCXX_THROW
#define _GLIBCXX_THROW(except) throw(except)
#endif

namespace rtt_memory
{
using namespace std;

unsigned total;
unsigned peak;

#if DRACO_DIAGNOSTICS & 2

bool is_active = false;

// We put the following in a wrapper so we can control destruction. We want to
// be sure is_active is forced to be false once alloc_map is destroyed.

struct memory_diagnostics
{
    map<void *, size_t> alloc_map;

    ~memory_diagnostics() { is_active = false; }
}
    st;

#endif // DRACO_DIAGNOSTICS & 2

//---------------------------------------------------------------------------------------//
bool set_memory_checking(bool new_status)
{
#if DRACO_DIAGNOSTICS & 2
    bool Result = is_active;

    total = 0;
    peak = 0;
    is_active = false;
    st.alloc_map.clear();
    is_active = new_status;
    
    return Result;
#endif
}

//---------------------------------------------------------------------------------------//
unsigned total_allocation()
{
    return total;
}

//---------------------------------------------------------------------------------------//
unsigned peak_allocation()
{
    return peak;
}

} // end namespace rtt_memory

using namespace rtt_memory;

#if DRACO_DIAGNOSTICS & 2
//---------------------------------------------------------------------------------------//
void *operator new(size_t n) _GLIBCXX_THROW(std::bad_alloc)
{
    void *Result = malloc(n);
    if (is_active)
    {
        total += n;
        // Don't use max() here; doing it with if statement allows programmers
        // to set a breakpoint here to find high water marks of memory usage.
        if (total>peak)
        {
            peak = total;
        }
        is_active = false;
        st.alloc_map[Result] = n;
        is_active = true;
    }
    return Result;
}

//---------------------------------------------------------------------------------------//
void operator delete(void *ptr) throw()
{
    free(ptr);
    if (is_active)
    {
        map<void *, size_t>::iterator i = st.alloc_map.find(ptr);
        if (i != st.alloc_map.end())
        {
            total -= i->second;
            is_active = false;
            st.alloc_map.erase(i);
            is_active = true;
        }
    }
}
#endif

//---------------------------------------------------------------------------//
// end of memory.cc
//---------------------------------------------------------------------------//
