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
 
void Timer::start()
{
    // Require(timerState == OFF);
    time -= mysec();
    timerState = ON;
    ++count;
}

void Timer::stop()
{
    // Require(timerState == ON);
    time += mysec();
    timerState = OFF;
}

void Timer::reset()
{
    Require(timerState == OFF);
    time = 0.0;
    count = 0;
}

Timer &Timer::operator+=(const Timer &rhs)
{
    Require(timerState == OFF);
    Require(rhs.timerState == OFF);
    time += rhs.time;
    count += rhs.count;
    return *this;
}

} // end namespace rtt_stopwatch

//---------------------------------------------------------------------------//
//                              end of Timer.cc
//---------------------------------------------------------------------------//
