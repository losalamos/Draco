//----------------------------------*-C++-*----------------------------------//
// Timer.cc
// Geoffrey Furnish
// Fri Oct  6 08:12:58 1995
//---------------------------------------------------------------------------//
// @> Stopwatch facility for GTS.
//---------------------------------------------------------------------------//

#include <iostream.h>
#include <stdio.h>
#include <string.h>

#include "util/Timer.hh"

//---------------------------------------------------------------------------//
// Constructor to use when the time will be recovered by calling the
// seconds() method.
//---------------------------------------------------------------------------//

Timer::Timer()
    : s(NULL), print_at_end(0)
{
    start = clock();
}

//---------------------------------------------------------------------------//
// Constructor which will cause the timer to print out the time when it goes
// out of scope.  Pass in a string suitable for use as the format string for
// sprintf, with a %f field to take the time in seconds.
//---------------------------------------------------------------------------//

Timer::Timer( char *_s )
    : s(_s), print_at_end(1)
{
    start = clock();
}

//---------------------------------------------------------------------------//
// Print out the final time, if requested.
//---------------------------------------------------------------------------//

Timer::~Timer()
{
    if (print_at_end) {
	int l = strlen(s);
	char *buf = new char[ l+20 ];

	sprintf( buf, s, seconds() );
	cout << buf << flush;

	delete[] buf;
    }
}

//---------------------------------------------------------------------------//
// Return the seconds accumulated on the clock.
//---------------------------------------------------------------------------//

float Timer::seconds()
{
    return (clock()-start)/float(CLOCKS_PER_SEC);
}

//---------------------------------------------------------------------------//
//                              end of Timer.cc
//---------------------------------------------------------------------------//
