//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <tpkg>/<class>.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright © 2003 The Regents of the University of California.
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
	if (string(argv[arg]) == "--version")
	{
	    if (rtt_c4::node() == 0)
		cout << argv[0] << ": version " << rtt_<pkg>::release() 
		     << endl;
	    rtt_c4::finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing <class>, " << ass.what()
	     << endl;
	rtt_c4::finalize();
	return 1;
    }

    {
	rtt_c4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_<spkg>::passed) 
	{
	    cout << "**** <class> Test: PASSED on " 
		 << rtt_c4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing <class> on " << rtt_c4::node() << endl;
    
    rtt_c4::finalize();
}   

//---------------------------------------------------------------------------//
//                        end of <class>.cc
//---------------------------------------------------------------------------//
