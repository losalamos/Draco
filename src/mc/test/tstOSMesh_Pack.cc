//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstOSMesh_Pack.cc
 * \author Thomas M. Evans
 * \date   Thu Dec 20 16:37:05 2001
 * \brief  OS_Mesh::pack function test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "MC_Test.hh"
#include "../OS_Mesh.hh"
#include "../Layout.hh"
#include "../XYCoord_sys.hh"
#include "../XYZCoord_sys.hh"
#include "../OS_Builder.hh"
#include "../Math.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_mc::XYCoord_sys;
using rtt_mc::XYZCoord_sys;
using rtt_mc::Layout;
using rtt_mc::OS_Mesh;
using rtt_mc::OS_Builder;
using rtt_mc::global::soft_equiv;
using rtt_mc_test::Parser;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void pack_2D()
{
    // make a builder from parser input 
    SP<Parser> parser(new Parser());
    OS_Builder builder(parser);

    SP<OS_Mesh> m1 = builder.build_Mesh();
    
    if (m1->num_cells() != 6) ITFAILS;

    // pack the mesh replicated
    SP<OS_Mesh::Pack> pack_m1 = m1->pack();
    if (pack_m1->get_num_packed_cells() != 6) ITFAILS;
    
    // unpack the mesh
    SP<OS_Mesh> m2 = pack_m1->unpack();

    // check for equality 
    if (*m1 != *m2) ITFAILS;
}

//---------------------------------------------------------------------------//

void pack_3D()
{
    // make a builder from parser input 
    SP<Parser> parser(new Parser("OS_Input_3D"));
    OS_Builder builder(parser);

    SP<OS_Mesh> m1;
    SP<OS_Mesh> m2;
    SP<OS_Mesh> m3;

    // build the mesh
    m1 = builder.build_Mesh();
    if (m1->num_cells() != 12) ITFAILS;

    // pack the mesh replicated
    {
	SP<OS_Mesh::Pack> pack_m1 = m1->pack();
	if (pack_m1->get_num_packed_cells() != 12) ITFAILS;
    
	// unpack the mesh
	m2 = pack_m1->unpack();
    }

    // check for equality 
    if (*m1 != *m2) ITFAILS;

    // pack the mesh using one-to-one mapping
    {
	vector<int> cell_list(m1->num_cells());
	for (int i = 0; i < cell_list.size(); i++)
	    cell_list[i] = i+1;
	SP<OS_Mesh::Pack> pack_m1 = m1->pack(cell_list);
	if (pack_m1->get_num_packed_cells() != 12) ITFAILS;

	// unpack the mesh
	m3 = pack_m1->unpack();
    }

    // m1 and m3 should be equivalent but not equal
    if (*m1 == *m3) ITFAILS;
    {
	if (m1->num_cells() != m3->num_cells())         ITFAILS;

	// the coordinate systems should be the same
	if (m3->get_Coord().get_Coord() != "xyz")       ITFAILS;

	// the layouts should be equal
	if (m1->get_Layout() != m3->get_Layout())       ITFAILS;
    
	// the vertices and cell pairs will be equivalent, but different
	if (m1->get_vertex() == m3->get_vertex())       ITFAILS;
	if (m1->get_cell_pair() == m3->get_cell_pair()) ITFAILS;

	// check the vertices for duplication
	OS_Mesh::vf_double vertex = m3->get_vertex();
	if (vertex.size() != 3) ITFAILS;
	for (int i = 0; i < vertex.size(); i++)
	    if (vertex[i].size() != 36) ITFAILS;

	// check that the meshes are equivalent
	double m1_dim, m3_dim;
	double m1_pos, m3_pos;
	for (int d = 1; d <= 3; d++)
	    for (int cell = 1; cell <= m1->num_cells(); cell++)
	    {
		m1_dim = m1->dim(d, cell);
		m3_dim = m3->dim(d, cell);
		if (!soft_equiv(m1_dim, m3_dim)) ITFAILS;

		m1_pos = m1->pos(d, cell);
		m3_pos = m3->pos(d, cell);
		if (!soft_equiv(m1_pos, m3_pos)) ITFAILS;
	    }
    }

    if (!m1->full_Mesh()) ITFAILS;
    if (m2->full_Mesh())  ITFAILS;
    if (m3->full_Mesh())  ITFAILS;
}

//---------------------------------------------------------------------------//

void pack_3D_compress()
{   
    // make a builder from parser input 
    SP<Parser> parser(new Parser("OS_Input_3D"));
    OS_Builder builder(parser);

    SP<OS_Mesh> m1;
    SP<OS_Mesh> m2;

    m1 = builder.build_Mesh();

    // pack up the mesh from 12 cells to 3 cells (see pg 139, vol II)
    {
	vector<int> cell_list(m1->num_cells());
	cell_list[0]  = 1;
	cell_list[1]  = -1;
	cell_list[2]  = 0;
	cell_list[3]  = 2;
	cell_list[4]  = -2;
	cell_list[5]  = 0;
	cell_list[6]  = 3;
	cell_list[7]  = -3;
	cell_list[8]  = 0;
	cell_list[9]  = -4;
	cell_list[10] = 0;
	cell_list[11] = 0;

	SP<OS_Mesh::Pack> pack = m1->pack(cell_list);
	if (pack->get_num_packed_cells() != 3) ITFAILS;

	// unpack the new mesh
	m2 = pack->unpack();
    }

    // now check the mesh
    {
	if (m2->num_cells() != 3) ITFAILS;
	
	// check the vertices for duplication
	OS_Mesh::vf_double vertex = m2->get_vertex();
	if (vertex.size() != 3) ITFAILS;
	for (int d = 0; d < 3; d++)
	    if (vertex[d].size() != 16) ITFAILS;

	// check the positions and dimensions of each cell
	if (m2->pos(1,1) != -0.5) ITFAILS;
	if (m2->pos(2,1) != 0.0)  ITFAILS;
	if (m2->pos(3,1) != 0.5)  ITFAILS;

	if (m2->pos(1,2) != -0.5) ITFAILS;
	if (m2->pos(2,2) != 2.0)  ITFAILS;
	if (m2->pos(3,2) != 0.5)  ITFAILS;

	if (m2->pos(1,3) != -0.5) ITFAILS;
	if (m2->pos(2,3) != 0.0)  ITFAILS;
	if (m2->pos(3,3) != 1.5)  ITFAILS;

	if (m2->dim(1,1) != 1.0)  ITFAILS;
	if (m2->dim(2,1) != 2.0)  ITFAILS;
	if (m2->dim(3,1) != 1.0)  ITFAILS;

	if (m2->dim(1,2) != 1.0)  ITFAILS;
	if (m2->dim(2,2) != 2.0)  ITFAILS;
	if (m2->dim(3,2) != 1.0)  ITFAILS;

	if (m2->dim(1,3) != 1.0)  ITFAILS;
	if (m2->dim(2,3) != 2.0)  ITFAILS;
	if (m2->dim(3,3) != 1.0)  ITFAILS;

	if (m2->full_Mesh())      ITFAILS;
    }

    // try to pack this again with a different mapping and we should get an
    // assertion
    bool caught = false;
    {
	vector<int> cell_list(3);
	cell_list[0] = 1;
	cell_list[1] = 3;
	cell_list[2] = 2;

	try
	{
	    SP<OS_Mesh::Pack> pack = m2->pack(cell_list);
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

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // this is a serial test
    if (C4::node())
    {
	C4::Finalize();
	return 0;
    }

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

	pack_2D();
	pack_3D();
	pack_3D_compress();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstOSMesh_Pack, " << ass.what()
	     << endl;
	C4::Finalize();
	return 1;
    }

    {
	// status of test
	cout << endl;
	cout <<     "*********************************************" << endl;
	if (rtt_mc_test::passed) 
	{
	    cout << "**** tstOSMesh_Pack Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }

    cout << "Done testing tstOSMesh_Pack on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstOSMesh_Pack.cc
//---------------------------------------------------------------------------//
