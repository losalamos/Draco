//----------------------------------*-C++-*----------------------------------//
// Stopwatch.hh
// Shawn Pautz
// Tue Mar 23 16:12:23 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifdef SDP

#ifndef __stopwatch_Stopwatch_hh__
#define __stopwatch_Stopwatch_hh__

#include "../ds++/NestMap.hh"
#include "Timer.hh"
#include <string>

using rtt_ds::NestMap;
using std::string;

namespace rtt_stopwatch
{
 
//===========================================================================//
// class Stopwatch - 
//
// Purpose :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Stopwatch 
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    typedef NestMap<string, Timer>::iterator iterator;

  private:

    // DATA

    NestMap<string, Timer> timers;
    
  public:

    // CREATORS
    
    Stopwatch() : timers() {}
    ~Stopwatch() {}

    // MANIPULATORS

    void start(const string& key);
    void stop();
    void stop(const string& key);
    
    // ACCESSORS

    iterator begin() { return timers.begin(); }
    iterator end() { return timers.end(); }

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_stopwatch

#endif                          // __stopwatch_Stopwatch_hh__

#endif

//---------------------------------------------------------------------------//
//                              end of stopwatch/Stopwatch.hh
//---------------------------------------------------------------------------//
