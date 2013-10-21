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

// We put the following in a wrapper so we can control destruction. We want to
// be sure is_active is set to false once alloc_map is destroyed.

#if DRACO_DIAGNOSTICS & 2
struct memory_diagnostics
{
    unsigned total;
    unsigned peak;
    bool is_active = false;
    map<void *, size_t> alloc_map;

    ~memory_diagnostics() { is_active = false; }
}
    st;

//---------------------------------------------------------------------------------------//
bool set_memory_checking(bool new_status)
{
    bool Result = st.is_active;

    st.total = 0;
    st.peak = 0;
    st.is_active = false;
    st.alloc_map.clear();
    st.is_active = new_status;
    
    return Result;
}

//---------------------------------------------------------------------------------------//
unsigned total_allocation()
{
    return st.total;
}

//---------------------------------------------------------------------------------------//
unsigned peak_allocation()
{
    return st.peak;
}
#endif

} // end namespace rtt_memory

using namespace rtt_memory;

#if DRACO_DIAGNOSTICS & 2
//---------------------------------------------------------------------------------------//
void *operator new(size_t n) _GLIBCXX_THROW(std::bad_alloc)
{
    void *Result = malloc(n);
    if (st.is_active)
    {
        st.total += n;
        st.peak = max(st.peak, st.total);
        st.is_active = false;
        st.alloc_map[Result] = n;
        st.is_active = true;
    }
    return Result;
}

//---------------------------------------------------------------------------------------//
void operator delete(void *ptr) throw()
{
    free(ptr);
    if (st.is_active)
    {
        map<void *, size_t>::iterator i = st.alloc_map.find(ptr);
        if (i != st.alloc_map.end())
        {
            st.total -= i->second;
            st.is_active = false;
            st.alloc_map.erase(i);
            st.is_active = true;
        }
    }
}
#endif

//---------------------------------------------------------------------------//
// end of memory.cc
//---------------------------------------------------------------------------//
