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

using rtt_stopwatch::Timer;
using std::cout;
using std::endl;

bool passed = true;

void pause( int seconds )
{
    double time = clock()/CLOCKS_PER_SEC;
    while (clock()/CLOCKS_PER_SEC - time < seconds) {}
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

    Timer testTimer1;
    Timer testTimer2;
    int pauseLength = 5;

    testTimer1.start();
    pause(pauseLength);
    testTimer1.stop();
    if (fabs(testTimer1.getTime() - pauseLength) > 0.01)
        passed = false;
    if (testTimer1.getCount() != 1)
        passed = false;

    testTimer1.start();
    pause(pauseLength);
    testTimer2.start();
    pause(pauseLength);
    testTimer2.stop();
    pause(pauseLength);
    testTimer1.stop();
    if (fabs(testTimer1.getTime() - 4*pauseLength) > 0.01)
        passed = false;
    if (fabs(testTimer2.getTime() - pauseLength) > 0.01)
        passed = false;
    if (testTimer1.getCount() != 2)
        passed = false;
    if (testTimer2.getCount() != 1)
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
