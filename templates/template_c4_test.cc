//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <tpkg>/<class>.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "../Release.hh"
#include "<spkg>.hh"

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//



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
		     << <namespace>::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
    }
    catch (std::exception &err)
    {
	std::cout << "ERROR: While testing <class>, " 
		  << err.what()
		  << std::endl;
	rtt_c4::finalize();
	return 1;
    }
    catch( ... )
    {
	std::cout << "ERROR: While testing <class>, " 
		  << "An unknown exception was thrown on processor "
                  << rtt_c4::node() << std::endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	std::cout << std::endl;
	std::cout <<     "*********************************************" 
                  << std::endl;
	if (<namespace>_test::passed) 
	{
	    std::cout << "**** <class> Test: PASSED on " 
                      << rtt_c4::node() 
                      << std::endl;
	}
	std::cout <<     "*********************************************" 
                  << std::endl;
	std::cout << std::endl;
    }
    
    rtt_c4::global_barrier();

    std::cout << "Done testing <class> on " << rtt_c4::node() 
              << std::endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of <class>.cc
//---------------------------------------------------------------------------//
