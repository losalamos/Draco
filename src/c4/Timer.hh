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
// 
//===========================================================================//

class Timer 
{
  private:
    // Beginning wall clock time.
    clock_t begin;

    // Flag determining if begin has been set.
    bool begin_set;
    
    // Ending wall clock time.
    clock_t end;

    // POSIX tms structure for beginning time.
    tms tms_begin;

    // POSIX tms structure for ending time.
    tms tms_end;

    // System clock resolution.
    const long clock_tick;

  public:
    //! Constructor.
    Timer() : begin_set(false), clock_tick(sysconf(_SC_CLK_TCK)) {/*...*/}

    //! Set the beginning of time cycle.
    inline void start();

    //! Set the end of time cycle.
    void stop() { end = times(&tms_end); }

    //! Do lap times without resetting the timer.
    void lap_start() { if (!begin_set) start(); }

    //! Return the wall clock time in seconds.
    inline double wall_clock() const;

    //! Return the system cpu time in seconds.
    inline double system_cpu() const;

    //! Return the user cpu time in seconds.
    inline double user_cpu() const;

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
    begin     = times(&tms_begin);
    begin_set = true; 
}

//---------------------------------------------------------------------------//

double Timer::wall_clock() const
{
    return (end - begin) / static_cast<double>(clock_tick);
}

//---------------------------------------------------------------------------//

double Timer::system_cpu() const
{
    return (tms_end.tms_stime - tms_begin.tms_stime) /
	static_cast<double>(clock_tick);
}

//---------------------------------------------------------------------------//

double Timer::user_cpu() const
{
    return (tms_end.tms_utime - tms_begin.tms_utime) /
	static_cast<double>(clock_tick); 
}

} // end namespace rtt_c4

#endif                          // __c4_Timer_hh__

//---------------------------------------------------------------------------//
//                              end of c4/Timer.hh
//---------------------------------------------------------------------------//
