//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4_test/test/tstAbort.cc
 * \author Thomas M. Evans
 * \date   Thu Jun  2 09:28:02 2005
 * \brief  C4 Abort test.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "../global.hh"
#include "../SpinLock.hh"
#include "../Release.hh"
#include "c4_test.hh"

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void abort_test()
{
    cout << "Entering abort on proc " << rtt_c4::node() << endl;
    
    rtt_c4::global_barrier();

    // only abort from processor 0 for nice output
    if (rtt_c4::node() == 0)
    {
        cout << "Aborting from processor 0" << endl;
        rtt_c4::abort();
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (std::string(argv[arg]) == "--version")
	{
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " 
		     << rtt_c4::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	std::cout << std::endl;
	std::cout <<     "*********************************************" 
                  << std::endl;
	if (rtt_c4_test::passed) 
	{
	    std::cout << "**** tstAbort Test: PASSED on " 
                      << rtt_c4::node() 
                      << std::endl;
	}
	std::cout <<     "*********************************************" 
                  << std::endl;
	std::cout << std::endl;
    }
    
    rtt_c4::global_barrier();

    // run test here so we get a pass message; this should simply abort the
    // program at this point
    abort_test();

    std::cout << "Done testing tstAbort on " << rtt_c4::node() 
              << std::endl;
    
    rtt_c4::finalize();

    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstAbort.cc
//---------------------------------------------------------------------------//
