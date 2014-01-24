//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4_test/test/tstAbort.cc
 * \author Thomas M. Evans
 * \date   Thu Jun  2 09:28:02 2005
 * \brief  C4 Abort test.
 * \note   Copyright (C) 2005-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../global.hh"
#include "../SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include <iostream>

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
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::initialize(argc, argv);

    bool runtest = false;

    // version tag
    for (int arg = 1; arg < argc; arg++)
    {
	if (std::string(argv[arg]) == "--version")
	{
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " 
		     << rtt_dsxx::release()
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}
        else if (std::string(argv[arg]) == "--runtest")
        {
            runtest = true;
        }
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	std::cout << std::endl;
	std::cout <<     "*********************************************" 
                  << std::endl;
        std::cout << "**** tstAbort Test: PASSED on " 
                  << rtt_c4::node() 
                  << std::endl;
	std::cout <<     "*********************************************" 
                  << std::endl;
	std::cout << std::endl;
    }
    
    rtt_c4::global_barrier();

    // run test here so we get a pass message; this should simply abort the
    // program at this point;  only run the test if --runtest is given so we
    // don't hose the automated testing
    if (runtest)
        abort_test();

    std::cout << "Done testing tstAbort on " << rtt_c4::node() 
              << std::endl;
    
    rtt_c4::finalize();

    return 0;
}   

//---------------------------------------------------------------------------//
// end of tstAbort.cc
//---------------------------------------------------------------------------//
