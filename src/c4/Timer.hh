//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Timer.hh
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:35:07 2002
 * \brief  Timer class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __c4_Timer_hh__
#define __c4_Timer_hh__

#include <iostream>
#include "ds++/Assert.hh"
#include "C4_Functions.hh"

namespace rtt_c4
{
 
//===========================================================================//
/*!
 * \class Timer
 *
 * \brief POSIX standard timer.
 *
 * The Timer class is used to calculate wall clock, user cpu, and system cpu
 * timings.  It uses the POSIX standard times function, so it should work
 * well on all (POSIX) systems.
 *
 * The POSIX implementation of timers are described in Sec. 8.15 of "Advanced
 * Programming in the UNIX Environment" by Stevens.
 *
 * Usage:
 * \code
 * #include <iostream>
 * #include "c4/Timer.hh"
 * using rtt_c4::Timer;
 *
 * Timer t;
 * t.start();
 * // do stuff
 * t.stop();
 * std::cout << t.wall_clock() << std::endl;
 * \endcode
 *
 * \example c4/test/tstTime.cc
 */
// revision history:
// -----------------
// 0) original
// 1) 2003/01/21 Added sum_* member functions (Lowrie).
// 
//===========================================================================//

class Timer 
{
  private:
    //! Beginning wall clock time.
    double begin;
    
    //! Ending wall clock time.
    double end;

    //! POSIX tms structure for beginning time.
    tms tms_begin;

    //! POSIX tms structure for ending time.
    tms tms_end;

    //! The number of clock ticks per second.
    //! \sa man times
    int const posix_clock_ticks_per_second;

    //! Flag determining if timer is currently on.
    bool timer_on;

    //! True if we can access MPI timers.
    bool const isMPIWtimeAvailable;
    
    //! sum of wall clock time over all intervals.
    double sum_wall;

    //! sum of system clock time over all intervals.
    double sum_system;

    //! sum of system clock time over all intervals.
    double sum_user;

    //! number of time intervals.
    int num_intervals;

    //! determine if MPI Wtime is available.
    bool setIsMPIWtimeAvailable() const;
    
  public:
    
    Timer();
    virtual ~Timer() { /* empty */ };
    Timer( Timer const & rhs );
    inline void start();
    inline void stop();
    inline double wall_clock() const;
    inline double system_cpu() const;
    inline double user_cpu()   const;
    inline double posix_err()  const;

    //! Return the wall clock time in seconds, summed over all intervals.
    double sum_wall_clock() const { Require(! timer_on); return sum_wall; }

    //! Return the system cpu time in seconds, summed over all intervals.
    double sum_system_cpu() const { Require(! timer_on); return sum_system; }

    //! Return the user cpu time in seconds, summed over all intervals.
    double sum_user_cpu() const { Require(! timer_on); return sum_user; }

    //! Return the number of time intervals used in the sums.
    int intervals() const { Require(! timer_on); return num_intervals; }

    inline void reset();
    void print( std::ostream &, int p = 2 ) const;
};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

//! Set the beginning of the time interval.
void Timer::start()
{
    Require(! timer_on);
    timer_on = true;
    ++num_intervals;
    // set both begin and tms_begin.
    begin    = wall_clock_time( tms_begin );
}

//---------------------------------------------------------------------------//
//! Set the end of the time interval.
void Timer::stop()
{
    Require( timer_on );
    // set both end and tms_end.
    end      = wall_clock_time( tms_end );
    timer_on = false; 

    sum_wall   += wall_clock();
    sum_system += system_cpu();
    sum_user   += user_cpu();
}

//---------------------------------------------------------------------------//
//! Return the wall clock time in seconds, for the last interval.
double Timer::wall_clock() const
{
    Require(! timer_on);
    return( end - begin );
}

//---------------------------------------------------------------------------//
//! Return the system cpu time in seconds, for the last interval.
double Timer::system_cpu() const
{
    Require(! timer_on);
    return( tms_end.tms_stime - tms_begin.tms_stime )
	/ static_cast<double>(posix_clock_ticks_per_second);
}

//---------------------------------------------------------------------------//
//! Return the user cpu time in seconds, for the last interval.
double Timer::user_cpu() const
{
    Require(! timer_on);
    return( tms_end.tms_utime - tms_begin.tms_utime )
	/ static_cast<double>(posix_clock_ticks_per_second); 
}

//---------------------------------------------------------------------------//
//! The error in the posix timings
double Timer::posix_err() const
{
    return 1.0/static_cast<double>(posix_clock_ticks_per_second);
}

//---------------------------------------------------------------------------//
//! Reset the interval sums.
void Timer::reset()
{
    Require(! timer_on);

    begin         = 0.0;
    end           = 0.0;
    timer_on      = false;
    sum_wall      = 0.0;
    sum_system    = 0.0;
    sum_user      = 0.0;
    num_intervals = 0;
    return;
}

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//

inline std::ostream& operator<<( std::ostream &out, const Timer &t )
{
    t.print( out, 2 );
    return out;
}

} // end namespace rtt_c4


#endif                          // __c4_Timer_hh__

//---------------------------------------------------------------------------//
//                              end of c4/Timer.hh
//---------------------------------------------------------------------------//
