//----------------------------------*-C++-*----------------------------------//
// Stopwatch.cc
// Shawn Pautz
// Tue Mar 23 16:12:23 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "Stopwatch.hh"

#include <iostream>
using std::cout;
using std::endl;

namespace rtt_stopwatch
{
 
void Stopwatch::start(const string& key)
{
    iterator iter = timers.open(key);
    iter->second.start();
}

void Stopwatch::stop()
{
    (*(timers.current())).second.stop();
    timers.close();
}

void Stopwatch::stop(const string& key)
{
    if ((*(timers.current())).first != key)
	cout << "Invalid stop request" << endl;
    else
	stop();
}
    
} // end namespace rtt_stopwatch

//---------------------------------------------------------------------------//
//                              end of Stopwatch.cc
//---------------------------------------------------------------------------//
