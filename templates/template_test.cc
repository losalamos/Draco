//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <tpkg>/<class>.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "<spkg>.hh"
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
	    cout << argv[0] << ": version " << rtt_<pkg>::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing <class>, " << ass.what() << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_<spkg>::passed) 
    {
        cout << "**** <class> Test: PASSED on" << C4::node() 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing <class>." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of <class>.cc
//---------------------------------------------------------------------------//
