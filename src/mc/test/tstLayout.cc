//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstLayout.cc
 * \author Thomas M. Evans
 * \date   Wed Jul 19 17:42:55 2000
 * \brief  Layout class tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "MC_Test.hh"
#include "../AMR_Layout.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <deque>
#include <list>
#include <string>

using namespace std;

using rtt_mc::AMR_Layout;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__, __FILE__);

//---------------------------------------------------------------------------//
// AMR_Layout TESTS
//---------------------------------------------------------------------------//
// in these tests build the following mesh layout
/*
                  0
            _________________
            |    |    |  |  |
	    |    |    |5 |6 |
   grad = 0 | 1  | 2  |__|__| 0
	    |    |    |  |  |  
	    |    |    |3 |4 |
	    |____|____|__|__|

                  0
*/
//---------------------------------------------------------------------------//

void test_AMR_Layout_1()
{
    AMR_Layout layout;
    
    // set the number of cells
    layout.set_size(3);
    if (layout.num_cells() != 3) ITFAILS;

    // set the faces
    for (int i = 1; i <= layout.num_cells(); i++)
	layout.set_size(i, 4);

    for (int i = 1; i <= layout.num_cells(); i++)
	if (layout.num_faces(i) != 4) ITFAILS;
}

//---------------------------------------------------------------------------//

void test_AMR_Layout_2()
{
    AMR_Layout layout(3, 4);
}

//---------------------------------------------------------------------------//
// MAIN
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);
    
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    if (!C4::node())
		cout << argv[0] << ": version " << rtt_mc::release() 
		     << endl;
	    C4::Finalize();
	    return 0;
	}

    try
    {
	test_AMR_Layout_1();
	test_AMR_Layout_2();
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	cout << "Dumb ass you screwed up; assertion: " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "**********************************" 
	 << endl;
    if (passed) 
    {
        cout << "**** Layout Self Test: PASSED on " 
	     << C4::node() << endl;
    }
    cout <<     "**********************************" 
	 << endl;
    cout << endl;

    cout << "Done testing Layout on node " << C4::node() 
	 << endl;

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tstLayout.cc
//---------------------------------------------------------------------------//
