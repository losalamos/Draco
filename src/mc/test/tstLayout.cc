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
#include "../Layout.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <deque>
#include <list>
#include <string>
#include <algorithm>

using namespace std;

using rtt_mc::Layout;
using rtt_mc::AMR_Layout;
using rtt_dsxx::SP;

bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__, __FILE__);

//---------------------------------------------------------------------------//
// Layout TESTS
//---------------------------------------------------------------------------//
// in these tests build the following mesh layout
/*
                  0
            _______________
            |    |    |   |
	    | 4  | 5  | 6 |
   grad = 0 |____|____|___| 0
            |    |    |   |
	    | 1  | 2  | 3 |
	    |____|____|___|
	    
                  0

   compress to this

                  0
            ____________
            |    |    |   
	    | 3  | 4  | -2 
   grad = 0 |____|____|_
            |    |    |  
	    | 1  | 2  | -1
	    |____|____|_
	    
                  0
		  
*/
//---------------------------------------------------------------------------//

void test_Layout()
{
    // build the layout
    Layout layout(6,4);
    if (layout.num_cells() != 6) ITFAILS;
    {
	layout(1,1) = 1;
	layout(1,2) = 2;
	layout(1,3) = 0;
	layout(1,4) = 4;

	layout(2,1) = 1;
	layout(2,2) = 3;
	layout(2,3) = 0;
	layout(2,4) = 5;

	layout(3,1) = 2;
	layout(3,2) = 0;
	layout(3,3) = 0;
	layout(3,4) = 6;

	layout(4,1) = 4;
	layout(4,2) = 5;
	layout(4,3) = 1;
	layout(4,4) = 0;

	layout(5,1) = 4;
	layout(5,2) = 6;
	layout(5,3) = 2;
	layout(5,4) = 0;

	layout(6,1) = 5;
	layout(6,2) = 0;
	layout(6,3) = 3;
	layout(6,4) = 0;
    }

    // check some cells
    if (layout(5,1) != 4) ITFAILS;
    if (layout(1,2) != 2) ITFAILS;
    if (layout(2,3) != 0) ITFAILS;
    if (layout(6,4) != 0) ITFAILS;

    for (int i = 1; i <= layout.num_cells(); i++)
	if (layout.num_faces(i) != 4) ITFAILS;

    // pack the full layout
    SP<Layout> unpacked_full_layout;
    SP<Layout::Pack> pack_layout;
    {
	vector<int> list(6);
	for (int i = 0; i < 6; i++)
	    list[i] = i+1;
	
	pack_layout = layout.pack(list);

	// check copy constructor
	if (pack_layout->get_size() != 31)            ITFAILS;
	if (pack_layout->get_num_packed_cells() != 6) ITFAILS;

	Layout::Pack test_pack = *pack_layout;
	SP<Layout> test_unpack = test_pack.unpack();

	if (*test_unpack != layout)       ITFAILS;

	// check pack accessors
	int  size = pack_layout->get_size();
	int *data = new int[size];

	copy(pack_layout->begin(), pack_layout->end(), data);

	Layout::Pack unpack(size, data);
	SP<Layout> unpack_layout = unpack.unpack();

	if (*unpack_layout != layout) ITFAILS;
    }
    
    // unpack and compare
    unpacked_full_layout = pack_layout->unpack();
	
    if (layout != *unpacked_full_layout) ITFAILS;

    // unpack the compressed layout
    SP<Layout> unpacked_com_layout;
    {
	vector<int> list(6);
	list[0] = 1;
	list[1] = 2;
	list[2] = -1;
	list[3] = 3;
	list[4] = 4;
	list[5] = -2;
	
	pack_layout = layout.pack(list);

	if (pack_layout->get_num_packed_cells() != 4) ITFAILS;
    }

    // unpack and check
    unpacked_com_layout = pack_layout->unpack();
    if (layout == *unpacked_com_layout) ITFAILS;

    // check the new layout
    Layout &lay = *unpacked_com_layout;
    if (lay.num_cells() != 4) ITFAILS;

    if (lay(1,1) != 1)  ITFAILS;
    if (lay(1,2) != 2)  ITFAILS;
    if (lay(1,3) != 0)  ITFAILS;
    if (lay(1,4) != 3)  ITFAILS;      

    if (lay(2,1) != 1)  ITFAILS;
    if (lay(2,2) != -1) ITFAILS;
    if (lay(2,3) != 0)  ITFAILS;
    if (lay(2,4) != 4)  ITFAILS;   

    if (lay(3,1) != 3)  ITFAILS;
    if (lay(3,2) != 4)  ITFAILS;
    if (lay(3,3) != 1)  ITFAILS;
    if (lay(3,4) != 0)  ITFAILS;  

    if (lay(4,1) != 3)  ITFAILS;
    if (lay(4,2) != -2) ITFAILS;
    if (lay(4,3) != 2)  ITFAILS;
    if (lay(4,4) != 0)  ITFAILS;

    // pack this layout and test
    SP<Layout> unpacked_again_com_layout;
    {
	vector<int> list(4);
	list[0] = 1;
	list[1] = 2;
	list[2] = 3;
	list[3] = 4;
	
	pack_layout = unpacked_com_layout->pack(list);
    }
    unpacked_again_com_layout = pack_layout->unpack();

    if (*unpacked_again_com_layout != *unpacked_com_layout) ITFAILS;

    // SPs should not be equal
    if (unpacked_again_com_layout == unpacked_com_layout)   ITFAILS;
}

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

   compress to this

                  0
            ____________
            |    |    |  
	    |    |    |-2
   grad = 0 | -1 | 1  |____
	    |    |    |  | 
	    |    |    |2 |-3
	    |____|____|__|_

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

    // pack the full layout
    SP<AMR_Layout::Pack> fullpack;
    SP<AMR_Layout>       full_unpacked;
    {
	vector<int> list(6);
	for (int i = 0; i < 6; i++)
	    list[i] = i+1;

	fullpack = layout.pack(list);
	
	if (fullpack->get_num_packed_cells() != 6)   ITFAILS;
	if (fullpack->get_size() != 1 + 6 + 24 + 25) ITFAILS;

	// check copy constructor
	AMR_Layout::Pack copy_pack = *fullpack;
	SP<AMR_Layout> copy_unpack = copy_pack.unpack();
	if (*copy_unpack != ref) ITFAILS;

	// check pack accessors
	int  size = fullpack->get_size();
	int *data = new int[size];
	copy(fullpack->begin(), fullpack->end(), data);

	AMR_Layout::Pack acc_pack(size, data);
	SP<AMR_Layout>   acc_unpack = acc_pack.unpack();
	if (*acc_unpack != ref) ITFAILS;
    }

    // unpack and compare
    full_unpacked = fullpack->unpack();
    if (*full_unpacked != ref) ITFAILS;

    // pack and unpack a compressed layout
    SP<AMR_Layout::Pack> compress_pack;
    SP<AMR_Layout>       compressed;
    {
	vector<int> list(6);
	list[0] = -1;
	list[1] = 1;
	list[2] = 2;
	list[3] = -3;
	list[4] = -2;
	list[5] = 0;

	compress_pack = layout.pack(list);

	if (compress_pack->get_num_packed_cells() != 2) ITFAILS;
    }

    // unpack and check
    compressed = compress_pack->unpack();
    AMR_Layout &refcom = *compressed;
    {
	if (refcom.num_cells() != 2) ITFAILS;

	if (refcom.num_cells_across(1,1) != 1) ITFAILS;
	if (refcom.num_cells_across(1,2) != 2) ITFAILS;
	if (refcom.num_cells_across(1,3) != 1) ITFAILS;
	if (refcom.num_cells_across(1,4) != 1) ITFAILS;
	if (refcom.num_cells_across(2,1) != 1) ITFAILS;
	if (refcom.num_cells_across(2,2) != 1) ITFAILS;
	if (refcom.num_cells_across(2,3) != 1) ITFAILS;
	if (refcom.num_cells_across(2,4) != 1) ITFAILS;

	if (refcom(1,1,1) != -1) ITFAILS;
	if (refcom(1,2,1) != 2)  ITFAILS;
	if (refcom(1,2,2) != -2) ITFAILS;
	if (refcom(1,3,1) != 0)  ITFAILS;
	if (refcom(1,4,1) != 0)  ITFAILS;

	if (refcom(2,1,1) != 1)  ITFAILS;
	if (refcom(2,2,1) != -3) ITFAILS;
	if (refcom(2,3,1) != 0)  ITFAILS;
	if (refcom(2,4,1) != -2) ITFAILS;
    }

    // check that you can't compress a compressed layout
    SP<AMR_Layout> compressed_again;
    {
	vector<int> list(2);
	list[0] = 1;
	list[1] = 2;

	compress_pack = compressed->pack(list);
    }
    compressed_again = compress_pack->unpack();

    if (*compressed_again != *compressed) ITFAILS;

    // now check the assertion
    bool caught = false;
    {
	vector<int> list(2);
	list[0] = 2;
	list[1] = 1;
	
	try
	{
	    compress_pack = compressed->pack(list);
	}
	catch (const rtt_dsxx::assertion &ass)
	{
	    cout << "Should catch this: " << ass.what() << endl;
	    caught = true;
	}
    }

    if (!caught) ITFAILS;
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
	// layout tests
	test_Layout();
	
	// AMR layout tests
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
