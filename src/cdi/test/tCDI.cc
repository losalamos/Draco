//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/tCDI.cc
 * \author Thomas M. Evans
 * \date   Tue Oct  9 15:52:01 2001
 * \brief  CDI test executable.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_test.hh"
#include "../Release.hh"
#include "ds++/Assert.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_cdi::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tCDI, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_cdi_test::passed) 
    {
        cout << "**** tCDI Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tCDI." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tCDI.cc
//---------------------------------------------------------------------------//
