//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstTime.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:19:16 2002
 * \brief  Test timing functions in C4.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "c4_test.hh"
#include "../Release.hh"
#include "../global.hh"
#include "../Timer.hh"
#include "../SpinLock.hh"
#include "ds++/Assert.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_c4::wall_clock_time;
using rtt_c4::wall_clock_resolution;
using rtt_c4::Timer;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void wall_clock_test()
{
    Timer t;

    double begin = rtt_c4::wall_clock_time();
    t.start();
    
    for (int i = 0; i < 200000000; i++)
    {
    }

    double end = rtt_c4::wall_clock_time();
    t.stop();
    
    cout << end-begin << endl;
    cout << t.wall_clock() << endl;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (C4::node() == 0)
		cout << argv[0] << ": version " << rtt_c4::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	
	wall_clock_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstTime, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	C4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_c4_test::passed) 
	{
	    cout << "**** tstTime Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstTime on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstTime.cc
//---------------------------------------------------------------------------//
