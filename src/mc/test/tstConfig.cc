//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstConfig.cc
 * \author Thomas M. Evans
 * \date   Mon Jan 14 11:20:29 2002
 * \brief  config.h tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"

#include <mc/config.h>

#include <sstream>
#include <typeinfo>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void check_eight_byte_int()
{
    int eight_byte = sizeof(EIGHT_BYTE_INT_TYPE);
    if (eight_byte == 8)
    {
	ostringstream message;
	message << "Eight byte int type set to " 
		<< typeid(EIGHT_BYTE_INT_TYPE).name();
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Eight byte int type of " 
		<< typeid(EIGHT_BYTE_INT_TYPE).name() << " NOT equal to "
		<< "eight bytes!";
	FAILMSG(message.str());
    }
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
		cout << argv[0] << ": version " << rtt_mc::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS

	// make sure the eight byte int type is properly set
	check_eight_byte_int();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstConfig, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	C4::HTSyncSpinLock slock;

	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_mc_test::passed) 
	{
	    cout << "**** tstConfig Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstConfig on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstConfig.cc
//---------------------------------------------------------------------------//
