//----------------------------------*-C++-*----------------------------------//
// MC_Test.hh
// Thomas M. Evans
// Wed Apr 21 19:33:36 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Some free functions that we will use throughout the MC test suite
//---------------------------------------------------------------------------//

#ifndef __mc_test_MC_Test_hh__
#define __mc_test_MC_Test_hh__

#include "../Layout.hh"
#include "../XYCoord_sys.hh"
#include "../OS_Mesh.hh"
#include "ds++/SP.hh"
#include <vector>
#include <iostream>

namespace rtt_mc_test
{

//===========================================================================//
// FAILURE LIMIT
//===========================================================================//

bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    return false;
}

//===========================================================================//
// 2D Objects
//===========================================================================//

//---------------------------------------------------------------------------//
// 2D Layout of a 3x3 mesh

rtt_mc::Layout build_2DLayout()
{
    using rtt_mc::Layout;

    // size of our 3x3 layout
    Layout layout = 6;
    
    // build the layout for our six cell problem
    for (int i = 1; i <= layout.num_cells(); i++)
    {
	// dimension cells
	layout.set_size(i, 4);

	// assign
	layout(i, 1) = i - 1;
	layout(i, 2) = i + 1;
	layout(i, 3) = i - 3;
	layout(i, 4) = i + 3;
    }

    // do layout boundaries
    for (int i = 1; i <= 3; i++)
	layout(i, 3) = 0;
    for (int i = 4; i <= 6; i++)
	layout(i, 4) = 0;
    layout(1, 1) = 1;
    layout(4, 1) = 4;
    layout(3, 2) = 0;
    layout(6, 2) = 0;

    return layout;
}    

//---------------------------------------------------------------------------//
// build 2D mesh (3x3)

dsxx::SP<rtt_mc::OS_Mesh> build_2DMesh()
{
    using rtt_mc::OS_Mesh;
    using rtt_mc::XYCoord_sys;
    using rtt_mc::Layout;
    using dsxx::SP;
    using namespace std;

    // build the Coord_sys
    SP<XYCoord_sys> coord;
    coord = new XYCoord_sys();

    // build the Layout
    Layout layout = build_2DLayout();

    // build the vertices
    vector<vector<double> > vert(2);
    for (int i = 0; i < 2; i++)
	vert[i].resize(12);
    vert[0][0] = -1;
    vert[1][0] = -1;
    for (int i = 1; i < 4; i ++)
    {
	vert[0][i] = vert[0][i-1] + 1;
	vert[1][i] = -1;
    }
    for (int i = 4; i < 8; i++)
    {
	vert[0][i] = vert[0][i-4];
	vert[1][i] = 1;
    }
    for (int i = 8; i < 12; i++)
    {
	vert[0][i] = vert[0][i-4];
	vert[1][i] = 3;
    }

    // build the cell pair
    vector<vector<int> > cp(6);
    for (int i = 0; i < 6; i++)
	cp[i].resize(4);

    for (int i = 1; i <= 3; i++)
	for (int j = 1; j <= 2; j++)
	{
	    int cell       = 1 + (i-1) + 3*(j-1);
	    int ref_vertex = 1 + (i-1) + 4*(j-1);
	    
	    cp[cell-1][0] = ref_vertex;
	    cp[cell-1][1] = ref_vertex + 1;
	    cp[cell-1][2] = ref_vertex + 4;
	    cp[cell-1][3] = ref_vertex + 5;
	}

    // now make the Mesh
    SP<OS_Mesh> mesh(new OS_Mesh(coord, layout, vert, cp));
    return mesh;
}

} // end namespace rtt_mc_test

#endif                          // __mc_test_MC_Test_hh__

//---------------------------------------------------------------------------//
//                              end of mc/test/MC_Test.hh
//---------------------------------------------------------------------------//
