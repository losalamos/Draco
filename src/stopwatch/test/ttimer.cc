//----------------------------------*-C++-*----------------------------------//
// ttimer.cc
// Shawn Pautz
// Thu Mar 25 13:40:43 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "../Timer.hh"
#include <iostream>
#include <cmath>

using rtt_stopwatch::Timer;
using rtt_stopwatch::OFF;
using rtt_stopwatch::ON;
using std::cout;
using std::endl;
using std::fabs;

bool passed = true;

void pause( int seconds )
{
    const long clock_tick(sysconf(_SC_CLK_TCK));
    clock_t begin;
    tms tmsNow;

    begin = times(&tmsNow);
    while ((times(&tmsNow)-begin)/static_cast<double>(clock_tick) < seconds) {}
}

void version(const std::string &progname)
{
    std::string version = "1.0.0";
    cout << progname << ": version " << version << endl;
}

int main( int argc, char *argv[] )
{
    for (int arg=1; arg < argc; arg++)
    {
	if (std::string(argv[arg]) == "--version")
	{
	    version(argv[0]);
	    return 0;
	}
    }
    
    cout << "Initiating test of the Timer class.\n";

    int pauseLength = 5;

    Timer timer;
    cout << "Pausing for " << pauseLength << " seconds" << endl;
    timer.start();
    pause(pauseLength);
    timer.stop();
    cout << "wallClock, systemCPU, userCPU: " << timer.wallClock()
	 << "  " << timer.systemCPU() << "  " << timer.userCPU() << endl;
    if (fabs(timer.wallClock() - pauseLength) > 0.1)
        passed = false;
    if (fabs((timer.systemCPU() + timer.userCPU())/timer.wallClock() - 1.)
	> 0.1)
        passed = false;
    if (timer.getCount() != 1)
        passed = false;

    timer.reset();
    if (timer.wallClock() != 0.)
        passed = false;
    if (timer.systemCPU() != 0.)
        passed = false;
    if (timer.userCPU() != 0.)
        passed = false;
    if (timer.getCount() != 0)
        passed = false;

    if (timer.getState() != OFF)
        passed = false;
    timer.start();
    if (timer.getState() != ON)
        passed = false;
    timer.stop();
    if (timer.getState() != OFF)
        passed = false;

    timer.reset();
    cout << "Pausing for " << 2*pauseLength << " seconds" << endl;
    timer.start();
    pause(pauseLength);
    timer.stop();
    timer.start();
    pause(pauseLength);
    timer.stop();
    cout << "wallClock, systemCPU, userCPU: " << timer.wallClock()
	 << "  " << timer.systemCPU() << "  " << timer.userCPU() << endl;
    if (fabs(timer.wallClock() - 2*pauseLength) > 0.1)
        passed = false;
    if (fabs((timer.systemCPU() + timer.userCPU())/timer.wallClock() - 1.)
	> 0.1)
        passed = false;
    if (timer.getCount() != 2)
        passed = false;

// Print the status of the test.

    cout << endl;
    cout <<     "*********************************" << endl;
    if (passed) 
    {
        cout << "**** Timer Self Test: PASSED ****" << endl;
    }
    else
    {
        cout << "**** Timer Self Test: FAILED ****" << endl;
    }
    cout <<     "*********************************" << endl;
    cout << endl;

    cout << "Done testing Timer class.\n";

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of ttimer.cc
//---------------------------------------------------------------------------//
