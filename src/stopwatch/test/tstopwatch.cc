//----------------------------------*-C++-*----------------------------------//
// tstopwatch.cc
// Shawn Pautz
// Tue Mar 30 12:24:54 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "../Stopwatch.hh"
#include <iostream>

#ifdef SDP

using rtt_stopwatch::Stopwatch;
using std::cout;
using std::endl;

bool passed = true;

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

#endif

int main( int argc, char *argv[] )
{

#ifdef SDP

    for (int arg=1; arg < argc; arg++)
    {
	if (std::string(argv[arg]) == "--version")
	{
	    version(argv[0]);
	    return 0;
	}
    }
    
    cout << "Initiating test of the Stopwatch class.\n";

    Stopwatch timers;

    timers.start("one");   // stack: 1
    timers.stop("one");    // stack:
    timers.start("one");   // stack: 1
    timers.stop();         // stack:
    timers.start("two");   // stack: 2
    timers.stop();         // stack:

    timers.start("one");          // stack: 1
    timers.start("one one");      // stack: 1 11
    timers.start("one one one");  // stack: 1 11 111
    timers.stop();                // stack: 1 11
    timers.stop();                // stack: 1
    timers.start("one two");      // stack: 1 12
    timers.stop();                // stack: 1
    timers.start("one three");    // stack: 1 13
    timers.stop();                // stack: 1
    timers.stop();                // stack:
    timers.start("three");        // stack: 3
    timers.start("three three");  // stack: 3 33
    timers.stop();                // stack: 3
    timers.stop();                // stack:

    cout << "\nDumping objects" << endl;
    for (Stopwatch::iterator iter = timers.begin();
	 iter != timers.end(); ++iter)
    {
	cout << (*iter).first << endl;
    }

// Print the status of the test.

    cout << endl;
    cout <<     "*************************************" << endl;
    if (passed) 
    {
        cout << "**** Stopwatch Self Test: PASSED ****" << endl;
    }
    else
    {
        cout << "**** Stopwatch Self Test: FAILED ****" << endl;
    }
    cout <<     "*************************************" << endl;
    cout << endl;

    cout << "Done testing Stopwatch class.\n";

#endif

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tstopwatch.cc
//---------------------------------------------------------------------------//
