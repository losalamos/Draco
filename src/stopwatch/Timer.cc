//----------------------------------*-C++-*----------------------------------//
// Timer.cc
// Shawn Pautz
// Thu Mar 25 11:15:07 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "Timer.hh"

#include "ds++/Assert.hh"

namespace rtt_stopwatch
{

Timer::Timer() : clockTick(sysconf(_SC_CLK_TCK)), elapsed(0),
		 timerState(OFF), count(0)
{
    tmsElapsed.tms_utime = 0;
    tmsElapsed.tms_stime = 0;
    tmsElapsed.tms_cutime = 0;
    tmsElapsed.tms_cstime = 0;
}

void Timer::start()
{
    Require(timerState == OFF);
    timerState = ON;
    ++count;
    begin = times(&tmsBegin);
}

void Timer::stop()
{
    Require(timerState == ON);
    timerState = OFF;
    end = times(&tmsEnd);
    elapsed += (end-begin);
    tmsElapsed.tms_utime += tmsEnd.tms_utime - tmsBegin.tms_utime;
    tmsElapsed.tms_stime += tmsEnd.tms_stime - tmsBegin.tms_stime;
    tmsElapsed.tms_cutime += tmsEnd.tms_cutime - tmsBegin.tms_cutime;
    tmsElapsed.tms_cstime += tmsEnd.tms_cstime - tmsBegin.tms_cstime;
}

void Timer::reset()
{
    Require(timerState == OFF);
    elapsed = 0;
    tmsElapsed.tms_utime = 0;
    tmsElapsed.tms_stime = 0;
    tmsElapsed.tms_cutime = 0;
    tmsElapsed.tms_cstime = 0;
    count = 0;
}

} // end namespace rtt_stopwatch

//---------------------------------------------------------------------------//
//                              end of Timer.cc
//---------------------------------------------------------------------------//
