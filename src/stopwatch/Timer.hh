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

#include <ctime>

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

#ifndef CLOCKS_PER_SEC
#ifdef CLK_TCK
#define CLOCKS_PER_SEC CLK_TCK  // Some systems still use this older form
#endif // CLK_TCK
#endif // CLOCKS_PER_SEC

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1000000  // For non-ANSI systems, like SunOS
#endif // CLOCKS_PER_SEC

enum state {OFF, ON};

inline double mysec()
{
    return (double)clock()/CLOCKS_PER_SEC;
}

class Timer 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    double time;
    state timerState;
    int count;
    
  public:

    // CREATORS
    
    Timer() : time(0.), timerState(OFF), count(0) {}
    ~Timer() {}

    // MANIPULATORS

    void start();
    void stop();
    void reset();

    Timer &operator+=(const Timer &rhs);
    
    // ACCESSORS

    double getTime() const { return time; }
    state getState() const { return timerState; }
    int getCount() const { return count; }
    static inline double resolution() { return 10.0/CLOCKS_PER_SEC; }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_stopwatch

#endif                          // __stopwatch_Timer_hh__

//---------------------------------------------------------------------------//
//                              end of stopwatch/Timer.hh
//---------------------------------------------------------------------------//
