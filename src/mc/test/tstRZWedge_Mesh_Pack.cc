//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstRZWedge_Mesh_Pack.cc
 * \author Thomas M. Evans
 * \date   Thu Dec 20 16:42:38 2001
 * \brief  RZWedge_Mesh::pack test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "MC_Test.hh"
#include "../RZWedge_Mesh.hh"
#include "../XYZCoord_sys.hh"
#include "../AMR_Layout.hh"
#include "../RZWedge_Builder.hh"
#include "../Release.hh"
#include "../Math.hh"
#include "viz/Ensight_Translator.hh"
#include "rng/Rnd_Control.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

using rtt_mc_test::Parser;
using rtt_mc_test::make_RZWedge_Mesh_AMR;
using rtt_mc::XYZCoord_sys;
using rtt_mc::AMR_Layout;
using rtt_dsxx::SP;
using rtt_mc::RZWedge_Mesh;
using rtt_mc::RZWedge_Builder;
using rtt_mc::global::soft_equiv;
using rtt_mc::global::pi;
using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_pack()
{
    // make a builder from parsing the RZWedge input
    SP<Parser> parser(new Parser("RZWedge_Input"));
    RZWedge_Builder builder(parser);

    SP<RZWedge_Mesh> m1;
    SP<RZWedge_Mesh> m2;
    SP<RZWedge_Mesh> m3;

    // build the mesh
    m1 = builder.build_Mesh();
    if (m1->num_cells() != 6) ITFAILS;

    // replicate the mesh
    {
	SP<RZWedge_Mesh::Pack> pack = m1->pack();
	if (pack->get_num_packed_cells() != 6) ITFAILS;
	m2 = pack->unpack();
    }

    if (m1 == m2)   ITFAILS;
    if (*m1 != *m2) ITFAILS;

    // replicate the mesh with a one-to-one pack
    {
	vector<int> cell_list(6);
	for (int i = 0; i < 6; i++)
	    cell_list[i] = i+1;

	SP<RZWedge_Mesh::Pack> pack = m1->pack(cell_list);
	if (pack->get_num_packed_cells() != 6) ITFAILS;
	m3 = pack->unpack();
    }
    
    if (m1 == m3)   ITFAILS;
    if (*m1 != *m3) ITFAILS;
}

//---------------------------------------------------------------------------//

void test_pack_AMR()
{
    SP<RZWedge_Mesh> m1;
    SP<RZWedge_Mesh> m2;
    SP<RZWedge_Mesh> m3;

    // make the amr mesh
    m1 = make_RZWedge_Mesh_AMR(5.0);
    
    // pack this replicated
    {
	SP<RZWedge_Mesh::Pack> pack = m1->pack();
	if (pack->get_num_packed_cells() != 12) ITFAILS;
	m2 = pack->unpack();
    }

    if (m1 == m2)   ITFAILS;
    if (*m1 != *m2) ITFAILS;

    // pack this compressed (save 5 cells, see pg 139 vol II)
    {
	vector<int> cell_list(12);
	cell_list[0]  = 1;
	cell_list[1]  = 2;
	cell_list[2]  = 3;
	cell_list[3]  = 4;
	cell_list[4]  = -1;
	cell_list[5]  = -3;
	cell_list[6]  = 5;
	cell_list[7]  = -2;
	cell_list[8]  = 0;
	cell_list[9]  = 0;
	cell_list[10] = 0;
	cell_list[11] = 0;

	SP<RZWedge_Mesh::Pack> pack = m1->pack(cell_list);
	if (pack->get_num_packed_cells() != 5) ITFAILS;
	m3 = pack->unpack();
    }
    
    // check compressed mesh
    if (*m1 == *m3)           ITFAILS;
    if (m3->num_cells() != 5) ITFAILS;
    {
	double pr    = 5.0 * pi / 180.0;
	double denom = sqrt(pr/sin(pr)) * cos(pr/2.0);

	// cell 1
	if (!soft_equiv(m3->get_low_x(1), 0.0))  ITFAILS;
	if (!soft_equiv(m3->get_high_x(1), 1.0)) ITFAILS;
	if (!soft_equiv(m3->get_low_z(1), 0.0))  ITFAILS;
	if (!soft_equiv(m3->get_high_z(1), 1.0)) ITFAILS;

	// cell 2
	if (!soft_equiv(m3->get_low_x(2), 0.0))  ITFAILS;
	if (!soft_equiv(m3->get_high_x(2), 0.5)) ITFAILS;
	if (!soft_equiv(m3->get_low_z(2), 1.0))  ITFAILS;
	if (!soft_equiv(m3->get_high_z(2), 1.5)) ITFAILS;

	// cell 3
	if (!soft_equiv(m3->get_low_x(3), 0.0))  ITFAILS;
	if (!soft_equiv(m3->get_high_x(3), 0.5)) ITFAILS;
	if (!soft_equiv(m3->get_low_z(3), 1.5))  ITFAILS;
	if (!soft_equiv(m3->get_high_z(3), 2.0)) ITFAILS;

	// cell 4
	if (!soft_equiv(m3->get_low_x(4), 0.5))  ITFAILS;
	if (!soft_equiv(m3->get_high_x(4), 1.0)) ITFAILS;
	if (!soft_equiv(m3->get_low_z(4), 1.0))  ITFAILS;
	if (!soft_equiv(m3->get_high_z(4), 1.5)) ITFAILS;

	// cell 5
	if (!soft_equiv(m3->get_low_x(5), 1.0))  ITFAILS;
	if (!soft_equiv(m3->get_high_x(5), 2.0)) ITFAILS;
	if (!soft_equiv(m3->get_low_z(5), 0.0))  ITFAILS;
	if (!soft_equiv(m3->get_high_z(5), 1.0)) ITFAILS;

	// do some next cell checks
	vector<double> r(3, 0.0);
	
	// cell 1
	if (m3->next_cell(1, 1) != 1)    ITFAILS;
	if (m3->next_cell(1, 2) != 5)    ITFAILS;
	if (m3->next_cell(1, 3) != 1)    ITFAILS;
	if (m3->next_cell(1, 4) != 1)    ITFAILS;
	if (m3->next_cell(1, 5) != 1)    ITFAILS;
	r[0] = 0.4;
	r[2] = 1.0;
	if (m3->next_cell(1, 6, r) != 2) ITFAILS;
	r[0] = 0.6;
	if (m3->next_cell(1, 6, r) != 4) ITFAILS;

	// cell 2
	if (m3->next_cell(2, 1) != 2)    ITFAILS;
	if (m3->next_cell(2, 2) != 4)    ITFAILS;
	if (m3->next_cell(2, 3) != 2)    ITFAILS;
	if (m3->next_cell(2, 4) != 2)    ITFAILS;
	if (m3->next_cell(2, 5) != 1)    ITFAILS;
	if (m3->next_cell(2, 6, r) != 3) ITFAILS; /* r has no effect,
						     only one cell across
						     face */

	// cell 3
	if (m3->next_cell(3, 1) != 3)    ITFAILS;
	if (m3->next_cell(3, 2) != -1)   ITFAILS;
	if (m3->next_cell(3, 3) != 3)    ITFAILS;
	if (m3->next_cell(3, 4) != 3)    ITFAILS;
	if (m3->next_cell(3, 5) != 2)    ITFAILS;
	if (m3->next_cell(3, 6) != -3)   ITFAILS;

	// cell 4
	if (m3->next_cell(4, 1) != 2)    ITFAILS;
	if (m3->next_cell(4, 2) != -2)   ITFAILS;
	if (m3->next_cell(4, 3) != 4)    ITFAILS;
	if (m3->next_cell(4, 4) != 4)    ITFAILS;
	if (m3->next_cell(4, 5) != 1)    ITFAILS;
	if (m3->next_cell(4, 6) != -1)   ITFAILS;

	// cell 5
	if (m3->next_cell(5, 1) != 1)    ITFAILS;
	if (m3->next_cell(5, 2) != 0)    ITFAILS;
	if (m3->next_cell(5, 3) != 5)    ITFAILS;
	if (m3->next_cell(5, 4) != 5)    ITFAILS;
	if (m3->next_cell(5, 5) != 5)    ITFAILS;
	if (m3->next_cell(5, 6) != -2)   ITFAILS;
    }

    if (!m1->full_Mesh()) ITFAILS;
    if (m2->full_Mesh())  ITFAILS;
    if (m3->full_Mesh())  ITFAILS;

    // try to pack again with a different mapping and we should get an
    // assertion
    bool caught = false;
    {
	vector<int> cell_list(5);
	cell_list[0] = 1;
	cell_list[1] = 2;
	cell_list[2] = 3;
	cell_list[3] = 5;
	cell_list[4] = 4;
	
	try
	{
	    SP<RZWedge_Mesh::Pack> pack = m3->pack(cell_list);
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

	test_pack();
	test_pack_AMR();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstRZWedge_Mesh_Pack, " << ass.what()
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
	    cout << "**** tstRZWedge_Mesh_Pack Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
    
    C4::gsync();

    cout << "Done testing tstRZWedge_Mesh_Pack on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstRZWedge_Mesh_Pack.cc
//---------------------------------------------------------------------------//
