//----------------------------------*-C++-*----------------------------------//
// Timer.hh
// Shawn Pautz
// Thu Mar 25 11:15:07 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __stopwatch_Timer_hh__
#define __stopwatch_Timer_hh__

#include <unistd.h>
#include <sys/times.h>

namespace rtt_stopwatch
{
 
//===========================================================================//
// class Timer - 
//
// Purpose : This is a basic timing class.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

enum state {OFF, ON};

class Timer 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    const long clockTick;

    clock_t begin;
    clock_t end;
    clock_t elapsed;

    tms tmsBegin;
    tms tmsEnd;
    tms tmsElapsed;

    state timerState;
    int count;

  public:

    // CREATORS

    Timer();

    // MANIPULATORS

    void start();
    void stop();
    void reset();

    double wallClock() const { return seconds(elapsed); }
    double systemCPU() const { return seconds(tmsElapsed.tms_stime); }
    double userCPU() const { return seconds(tmsElapsed.tms_utime); }
    state getState() const { return timerState; }
    int getCount() const { return count; }

  private:

    double seconds(const clock_t &time) const
    { return time/static_cast<double>(clockTick); }
    
    // IMPLEMENTATION
};

} // end namespace rtt_stopwatch

#endif                          // __stopwatch_Timer_hh__

//---------------------------------------------------------------------------//
//                              end of stopwatch/Timer.hh
//---------------------------------------------------------------------------//
