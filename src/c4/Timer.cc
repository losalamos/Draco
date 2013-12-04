//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Timer.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:56:11 2002
 * \note   Copyright (C) 2002-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "C4_sys_times.h"
#include "Timer.hh"
#include <iomanip>
#include <cstdlib>
#include <iostream>

namespace rtt_c4
{

#ifdef HAVE_PAPI

/* Initialize static non-const data members found in the Timer class */

// Count up the total L2 cache misses and hits and the total number of
// floating point operations. These allow us to report the percentage of cache
// hits and the number of floating point operations per cache miss.
int Timer::papi_events_[papi_max_counters_] = {
    PAPI_L2_TCM, PAPI_L2_TCH, PAPI_FP_OPS };

unsigned Timer::papi_num_counters_ = 0;

long long Timer::papi_raw_counts_[papi_max_counters_] = {
    0, 0, 0};

#endif

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//

//! Constructor
Timer::Timer()
    : begin( 0.0 ),
      end(   0.0 ),
      tms_begin( DRACO_TIME_TYPE() ),
      tms_end(   DRACO_TIME_TYPE() ),
      posix_clock_ticks_per_second( DRACO_CLOCKS_PER_SEC ),
      timer_on(  false ),
      isMPIWtimeAvailable( setIsMPIWtimeAvailable() ),
      sum_wall(      0.0 ),
      sum_system(    0.0 ),
      sum_user(      0.0 ),
      num_intervals( 0   )
{
#ifdef HAVE_PAPI
    
    // Initialize the PAPI library on construction of first timer.
    
    static bool first_time = true;
    if (first_time)
    {
	int retval = PAPI_library_init(PAPI_VER_CURRENT);
	if (retval != PAPI_VER_CURRENT)
        {
            std::cout << "PAPI library init error!" << std::endl;
            exit(EXIT_FAILURE);
	}
        
	if (PAPI_query_event(PAPI_FP_OPS) != PAPI_OK)
        {
	    std::cout << "PAPI: No floating operations counter" << std::endl;
        }

	if (PAPI_query_event(PAPI_L2_TCM) != PAPI_OK)
        {
	    std::cout << "PAPI: No L2 cache miss counter" << std::endl;
        }

	if (PAPI_query_event(PAPI_L2_TCH) != PAPI_OK)
        {
	    std::cout << "PAPI: No cache hit counter" << std::endl;
        }

	papi_num_counters_ = PAPI_num_counters();
        if (papi_num_counters_<sizeof(papi_events_)/sizeof(int))
        {
            std::cout << "PAPI: This system has only " << papi_num_counters_
                      << " hardware counters.\n" << std::endl;
            std::cout << "Some performance statistics will not be available."
                      << std::endl;
	}

        // At present, some platforms *lie* about how many counters they have
        // available, reporting they have three then returning an out of
        // counters error when you actually try to assign the three counter
        // types listed above. Until we have a fix, hardwire to leave out the
        // flops count, which is the least essential of the three counts.
        papi_num_counters_ = 2;

	if (papi_num_counters_ > sizeof(papi_events_)/sizeof(int))
            papi_num_counters_ = sizeof(papi_events_)/sizeof(int);

	int result = PAPI_start_counters(papi_events_, papi_num_counters_);
	if (result != PAPI_OK)
        {
	    std::cout << "Failed to start hardware counters with error "
                      << result << std::endl;
            exit(EXIT_FAILURE);
        }

        // sum of papi wall clock cycles
        sum_papi_wc_cycle = 0;
        // sum of papi wall clock time (microseconds)
        sum_papi_wc_usec = 0;
        // sum of papi virtual cycles
        sum_papi_virt_cycle = 0;
        // sum of papi virtual time (microseconds)
        sum_papi_virt_usec = 0;
        
	first_time = false;
    }
#endif

    reset();
}

//---------------------------------------------------------------------------//
// Member functions
//---------------------------------------------------------------------------//

//! Print out a timing report.
void Timer::print( std::ostream &out, int p ) const
{
    using std::setw;
    using std::ios;
    
    out.setf(ios::fixed, ios::floatfield);
    out.precision(p);
    out << '\n';
    
    if ( num_intervals > 1 )
	out << "LAST INTERVAL: " << '\n';

    out << setw(20) << "WALL CLOCK TIME: " << wall_clock() << " sec." << '\n';
    out << setw(20) << "  USER CPU TIME: " << user_cpu()   << " sec." << '\n';
    out << setw(20) << "SYSTEM CPU TIME: " << system_cpu() << " sec." << '\n';
    out << '\n';
    
    if ( num_intervals > 1 )
    {
	out << "OVER " << num_intervals << " INTERVALS: " << '\n';
	out << setw(20) << "WALL CLOCK TIME: " << sum_wall_clock()
	    << " sec." << '\n';
	out << setw(20) << "  USER CPU TIME: " << sum_user_cpu()
	    << " sec." << '\n';
	out << setw(20) << "SYSTEM CPU TIME: " << sum_system_cpu()
	    << " sec." << "\n\n";

#ifdef HAVE_PAPI
        double const miss = sum_L2_cache_misses();
        double const hit  = sum_L2_cache_hits();
        out << "PAPI Events:\n"

            << setw(26) << "L2 Cache misses  : "
            << sum_L2_cache_misses()     << "\n"

            << setw(26) << "L2 Cache hits    : "
            << sum_L2_cache_hits()       << "\n"

            << setw(26) << "Percent hit      : "
            << 100.0 * hit / (miss+hit)  << "\n"

            << setw(26) << "FP operations    : "
            << sum_floating_operations() << "\n"
            
            << setw(26) << "Wall Clock cycles: "
            << sum_papi_wc_cycles() << "\n"

            << setw(26) << "Wall Clock time (us): "
            << sum_papi_wc_usecs() << "\n"

            << setw(26) << "Virtual cycles: "
            << sum_papi_virt_cycles() << "\n"

            << setw(26) << "Virtual time (us): "
            << sum_papi_virt_usecs() << "\n"

            << std::endl;
#endif
    }
    
    out.flush();
}

//---------------------------------------------------------------------------//
//! Print out a timing report as a single line summary.
void Timer::printline( std::ostream &out,
                       unsigned const p,
                       unsigned const w ) const
{
    using std::setw;
    using std::ios;
    
    out.setf(ios::fixed, ios::floatfield);
    out.precision(p);

    // Width of first column (intervals) should be set by client before
    // calling this function.
    out << num_intervals
	<< setw(w) << sum_user_cpu()
	<< setw(w) << sum_system_cpu()
	<< setw(w) << sum_wall_clock();

#ifdef HAVE_PAPI
    double const miss = sum_L2_cache_misses();
    double const hit  = sum_L2_cache_hits();
    out << setw(w) 
	<< 100.0 * hit / (miss+hit);
    if (papi_num_counters_>2)
    {
	out << setw(w)
            << sum_floating_operations() / miss;
    }
#endif

    out << std::endl;

    out.flush();
}

//---------------------------------------------------------------------------//
// Is this an MPI or Posix timer?
//---------------------------------------------------------------------------//
bool Timer::setIsMPIWtimeAvailable() const
{
#ifdef C4_SCALAR
    return false;
#else
    return true;
#endif
}

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
// end of Timer.cc
//---------------------------------------------------------------------------//
