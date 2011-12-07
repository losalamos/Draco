//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Timer.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:56:11 2002
 * \brief  Timer member definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <c4/config.h>
#include "C4_sys_times.h"
#include "Timer.hh"
#include <iomanip>
#include <cstdlib>

namespace rtt_c4
{
  using std::cout;
  using std::endl;
  using std::min;

#ifdef HAVE_PAPI
/* static */
int Timer::papi_events_[papi_max_counters_] = {
    PAPI_L2_TCM, PAPI_L2_TCH, PAPI_FP_OPS };

/* static */
unsigned Timer::papi_num_counters_ = 0;

/* static */
long long Timer::papi_raw_counts_[papi_max_counters_] = {0, 0, 0};
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
      timer_on(          false ),
      isMPIWtimeAvailable( setIsMPIWtimeAvailable() ),
      sum_wall( 0.0 ),
      sum_system( 0.0 ),
      sum_user( 0.0 ),
      num_intervals( 0 )
{
#ifdef HAVE_PAPI
    // Initialize the PAPI library on construction of first timer
    static bool first_time = true;
    if (first_time)
    {
        int retval;
        // int EventSet = PAPI_NULL;
        // unsigned int native = 0x0;
        // PAPI_event_info_t info;
	
	/* Initialize the library */
	retval = PAPI_library_init(PAPI_VER_CURRENT);
	if (retval != PAPI_VER_CURRENT) {
            cout << "PAPI library init error!" << endl;
            exit(EXIT_FAILURE);
	}

	if (PAPI_query_event(PAPI_FP_OPS) != PAPI_OK)
        {
	    cout << "PAPI: No floating operations counter" << endl;
        }
	else
        {
	    cout << "PAPI: Floating operations are countable" << endl;
        }

	if (PAPI_query_event(PAPI_L2_TCM) != PAPI_OK)
        {
	    cout << "PAPI: No L2 cache miss counter" << endl;
        }
	else
        {
	    cout << "PAPI: L2 cache misses are countable" << endl;
        }

	if (PAPI_query_event(PAPI_L2_TCH) != PAPI_OK)
        {
	    cout << "PAPI: No cache hit counter" << endl;
        }
        else
        {
	    cout << "PAPI: L2 cache hits are countable" << endl;
        }

	papi_num_counters_ = PAPI_num_counters();
	cout << "PAPI: This system has " << papi_num_counters_
	     << " hardware counters." << endl;

	if (papi_num_counters_ > sizeof(papi_events_)/sizeof(int))
            papi_num_counters_ = sizeof(papi_events_)/sizeof(int);

	int result = PAPI_start_counters(papi_events_, papi_num_counters_);
	if (result != PAPI_OK)
        {
	    cout << "Failed to start hardware counteres with error "
		 << result << endl;
        }

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
    using std::endl;
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
	    << " sec." << '\n';
	out << '\n';
    }
    out.flush();
}

//! Print out a timing report as a single line summary.
  void Timer::printline( std::ostream &out,
			 unsigned const p,
			 unsigned const w ) const
{
    using std::setw;
    using std::endl;
    using std::ios;
    
    out.setf(ios::fixed, ios::floatfield);
    out.precision(p);

    // Width of first column (intervals) should be set by client
    // before calling this function.
    out << num_intervals
	<< setw(w) << sum_user_cpu()
	<< setw(w) << sum_system_cpu()
	<< setw(w) << sum_wall_clock();

#ifdef HAVE_PAPI
    double const miss = sum_L2_cache_misses();
    double const hit = sum_L2_cache_hits();
    out << setw(w) 
	<< 100.0 * hit / (miss+hit)
	<< setw(w)
	<< sum_floating_operations() / miss;
#endif

    out << endl;

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
//                              end of Timer.cc
//---------------------------------------------------------------------------//
