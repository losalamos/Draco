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

#include <unistd.h>
#include <sys/times.h>
#include <iostream>
#include "ds++/Assert.hh"

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
    // Beginning wall clock time.
    clock_t begin;
    
    // Ending wall clock time.
    clock_t end;

    // POSIX tms structure for beginning time.
    tms tms_begin;

    // POSIX tms structure for ending time.
    tms tms_end;

    // System clock resolution.
    const long clock_tick;

    // Flag determining if timer is currently on.
    bool timer_on;

    // sum of wall clock time between starts and stops
    double sum_wall;

    // sum of system clock time between starts and stops
    double sum_system;

    // sum of system clock time between starts and stops
    double sum_user;

  public:
    //! Constructor.
    Timer() : timer_on(false), clock_tick(sysconf(_SC_CLK_TCK)) { reset(); }

    //! Set the beginning of time cycle.
    inline void start();

    //! Set the end of time cycle.
    inline void stop();

    //! Return the wall clock time in seconds.
    inline double wall_clock() const;

    //! Return the system cpu time in seconds.
    inline double system_cpu() const;

    //! Return the user cpu time in seconds.
    inline double user_cpu() const;

    //! Return the wall clock time in seconds, summed over starts and stops.
    double sum_wall_clock() const { Require(! timer_on); return sum_wall; }

    //! Return the system cpu time in seconds, summed over starts and stops.
    double sum_system_cpu() const { Require(! timer_on); return sum_system; }

    //! Return the user cpu time in seconds, summed over starts and stops.
    double sum_user_cpu() const { Require(! timer_on); return sum_user; }

    //! Reset the sums of times.
    inline void reset();

    //! Print out a timing report.
    void print(std::ostream &, int p = 2) const;
};

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//

inline std::ostream& operator<<(std::ostream &out, const Timer &t)
{
    t.print(out, 2);
    return out;
}

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

void Timer::start()
{
    Require(! timer_on);
    timer_on = true; 
    begin    = times(&tms_begin);
}

//---------------------------------------------------------------------------//

void Timer::stop()
{
    Require(timer_on);
    timer_on = false; 
    end      = times(&tms_end);

    sum_wall   += wall_clock();
    sum_system += system_cpu();
    sum_user   += user_cpu();
}

//---------------------------------------------------------------------------//

double Timer::wall_clock() const
{
    Require(! timer_on);
    return (end - begin) / static_cast<double>(clock_tick);
}

//---------------------------------------------------------------------------//

double Timer::system_cpu() const
{
    Require(! timer_on);
    return (tms_end.tms_stime - tms_begin.tms_stime) /
	static_cast<double>(clock_tick);
}

//---------------------------------------------------------------------------//

double Timer::user_cpu() const
{
    Require(! timer_on);
    return (tms_end.tms_utime - tms_begin.tms_utime) /
	static_cast<double>(clock_tick); 
}

//---------------------------------------------------------------------------//

void Timer::reset()
{
    Require(! timer_on);

    sum_wall   = 0.0;
    sum_system = 0.0;
    sum_user   = 0.0;
}

} // end namespace rtt_c4

#endif                          // __c4_Timer_hh__

//---------------------------------------------------------------------------//
//                              end of c4/Timer.hh
//---------------------------------------------------------------------------//
