//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/test/tGandolfWithCDI.cc
 * \author Thomas M. Evans
 * \date   Mon Oct 29 17:16:32 2001
 * \brief  Gandolf + CDI test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_gandolf_test.hh"
#include "../Release.hh"
#include "ds++/Assert.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_gandolf_CDI()
{
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_cdi_gandolf::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	test_gandolf_CDI();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tGandolfWithCDI, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_cdi_gandolf_test::passed) 
    {
        cout << "**** tGandolfWithCDI Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tGandolfWithCDI." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tGandolfWithCDI.cc
//---------------------------------------------------------------------------//
