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

AMR_Layout test_AMR_Layout_1()
{
    AMR_Layout layout;
    
    // set the number of cells
    layout.set_size(6);
    if (layout.num_cells() != 6) ITFAILS;

    // set the faces
    for (int i = 1; i <= layout.num_cells(); i++)
	layout.set_size(i, 4);

    for (int i = 1; i <= layout.num_cells(); i++)
    {
	if (layout.num_faces(i) != 4) ITFAILS;

	for (int f = 1; f <= layout.num_faces(i); f++)
	    if (layout.num_cells_across(i,f) != 1) ITFAILS;
    }

    // build the layout
    layout(1,1,1) = 1;
    layout(1,2,1) = 2;
    layout(1,3,1) = 0;
    layout(1,4,1) = 0;

    layout.set_size(2,2,2);
    if (layout.num_cells_across(2,2) != 2) ITFAILS;

    layout(2,1,1) = 1;
    layout(2,2,1) = 3;
    layout(2,2,2) = 5;
    layout(2,3,1) = 0;
    layout(2,4,1) = 0;

    layout(3,1,1) = 2;
    layout(3,2,1) = 4;
    layout(3,3,1) = 0;
    layout(3,4,1) = 5;

    layout(4,1,1) = 3;
    layout(4,2,1) = 0;
    layout(4,3,1) = 0;
    layout(4,4,1) = 6;

    layout(5,1,1) = 2;
    layout(5,2,1) = 6;
    layout(5,3,1) = 3;
    layout(5,4,1) = 0;

    layout(6,1,1) = 5;
    layout(6,2,1) = 0;
    layout(6,3,1) = 4;
    layout(6,4,1) = 0;

    // check the layout
    vector<int> cells;
    {
	cells.resize(1);
	cells[0] = 2;
	if (layout(1,2) != cells) ITFAILS;

	cells[0] = 1;
	if (layout(1,1) != cells) ITFAILS;

	cells.resize(2);
	cells[0] = 3;
	cells[1] = 5;
	if (layout(2,2) != cells) ITFAILS;

	cells.resize(1);
	cells[0] = 6;
	if (layout(4,4) != cells) ITFAILS;
    }

    return layout;
}

//---------------------------------------------------------------------------//

void test_AMR_Layout_2(const AMR_Layout &ref)
{
    AMR_Layout layout(ref.num_cells(), 4);
    layout.set_size(2,2,2);

    layout(1,1,1) = 1;
    layout(1,2,1) = 2;
    layout(1,3,1) = 0;
    layout(1,4,1) = 0;

    layout(2,1,1) = 1;
    layout(2,2,1) = 3;
    layout(2,2,2) = 5;
    layout(2,3,1) = 0;
    layout(2,4,1) = 0;

    layout(3,1,1) = 2;
    layout(3,2,1) = 4;
    layout(3,3,1) = 0;
    layout(3,4,1) = 5;

    layout(4,1,1) = 3;
    layout(4,2,1) = 0;
    layout(4,3,1) = 0;
    layout(4,4,1) = 6;

    layout(5,1,1) = 2;
    layout(5,2,1) = 6;
    layout(5,3,1) = 3;
    layout(5,4,1) = 0;

    layout(6,1,1) = 5;
    layout(6,2,1) = 0;
    layout(6,3,1) = 4;
    layout(6,4,1) = 0;

    if (layout != ref) ITFAILS;
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
	AMR_Layout reference = test_AMR_Layout_1();
	test_AMR_Layout_2(reference);
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
