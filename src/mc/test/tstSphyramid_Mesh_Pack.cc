//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/tstSphyramid_Mesh_Pack.cc
 * \author Jeffery Densmore
 * \date   Mon Nov 24 11:45:03 2003
 * \brief  Sphyramid_Mesh::pack test.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc_test.hh"
#include "MC_Test.hh"
#include "../Sphyramid_Mesh.hh"
#include "../Sphyramid_Builder.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void test_pack()
{
    using rtt_mc_test::Parser;
    using rtt_dsxx::SP;
    using rtt_mc::Sphyramid_Mesh;
    using rtt_mc::Sphyramid_Builder;
    using std::vector;
    using rtt_mc::global::soft_equiv;
    using rtt_mc::global::pi;
    using std::cout;
    using std::endl;

    // make a builder from parsing the Sphyramid_Input file
    SP<Parser> parser(new Parser("Sphyramid_Input"));
    Sphyramid_Builder builder(parser);

    SP<Sphyramid_Mesh> m1;
    SP<Sphyramid_Mesh> m2;
    SP<Sphyramid_Mesh> m3;
    SP<Sphyramid_Mesh> m4;

    // build the mesh
    m1 = builder.build_Mesh();
    if (m1->num_cells() != 5) ITFAILS;

    // replicate the mesh
    {
	SP<Sphyramid_Mesh::Pack> pack = m1->pack();
	if (pack->get_num_packed_cells() != 5) ITFAILS;
	m2 = pack->unpack();
    }

    if ( m1 ==  m2) ITFAILS;
    if (*m1 != *m2) ITFAILS;

    // replicate the mesh with a one-to-one pack
    {
	vector<int> cell_list(5);
	for (int i = 0; i < 5; i++)
	{
	    cell_list[i] = i+1;
	}
	SP<Sphyramid_Mesh::Pack> pack = m1->pack(cell_list);
	if (pack->get_num_packed_cells() != 5) ITFAILS;
	m3 = pack->unpack();
    }

    if ( m1 ==  m3) ITFAILS;
    if (*m1 != *m3) ITFAILS;

    // replicate using only first and last cell (and reverse order)
    { 
	vector<int> cell_list(5);    
	cell_list[0] =  2;
	cell_list[1] = -1;
	cell_list[2] = -2;
	cell_list[3] = -3;
	cell_list[4] =  1;
    
	SP<Sphyramid_Mesh::Pack> pack = m1->pack(cell_list);
	if(pack->get_num_packed_cells() != 2) ITFAILS;
	m4 = pack->unpack();
    }

    // check compressed mesh
    if (*m1             == *m4) ITFAILS;
    if (m4->num_cells() != 2)   ITFAILS;
    {
	// spherical cone angle (alpha) is 45 degrees
	double r_to_x   = pow((2.*(1.-cos(pi/4.)))/(tan(pi/4)*tan(pi/4)), 1./3.);
       	
	// the radial boundaries are 0,0.5,1,9/7,13/7,3
	
	// cell 2 (was cell 1)
	if (!soft_equiv(m4->get_low_x(2) , 0.0))        ITFAILS;
	if (!soft_equiv(m4->get_high_x(2), 0.5*r_to_x)) ITFAILS;

	// cell 1 (was cell 5)
	if (!soft_equiv(m4->get_low_x(1) , 13./7.*r_to_x)) ITFAILS;
	if (!soft_equiv(m4->get_high_x(1), 3.0*r_to_x))    ITFAILS;

	// do some next cell checks
	
	// cell 1
	if(m4->next_cell(1,1) != -3)  ITFAILS;
	if(m4->next_cell(1,2) !=  0)  ITFAILS;
	if(m4->next_cell(1,3) !=  1)  ITFAILS;
	if(m4->next_cell(1,4) !=  1)  ITFAILS;
	if(m4->next_cell(1,5) !=  1)  ITFAILS;
	if(m4->next_cell(1,6) !=  1)  ITFAILS;
	
	// cell 2
	if(m4->next_cell(2,1) !=  2)  ITFAILS;
	if(m4->next_cell(2,2) != -1)  ITFAILS;
	if(m4->next_cell(2,3) !=  2)  ITFAILS;
	if(m4->next_cell(2,4) !=  2)  ITFAILS;
	if(m4->next_cell(2,5) !=  2)  ITFAILS;
	if(m4->next_cell(2,6) !=  2)  ITFAILS;

	
    }
    
    if (!m1->full_Mesh()) ITFAILS;
    if (m2->full_Mesh())  ITFAILS;
    if (m3->full_Mesh())  ITFAILS;
    if (m4->full_Mesh())  ITFAILS;

    // try to pack again with a different mapping and we should get an 
    // assertion
    bool caught = false;
    {
	vector<int> cell_list(2);
	cell_list[0] = 2;
	cell_list[1] = 1;

	try
	{
	    SP<Sphyramid_Mesh::Pack> pack = m4->pack(cell_list);
	}
	catch (const rtt_dsxx::assertion &ass)
	{
	    cout << "Should catch this: " << ass.what() << endl;
	    caught = true;
	}
    }
    if(!caught) ITFAILS;
    
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    using std::string;
    using std::cout;
    using std::endl;
    
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


    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstSphyramid_Mesh_Pack.cc, " << ass.what()
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
	    cout << "**** tstSphyramid_Mesh_Pack.cc Test: PASSED on " 
		 << C4::node() << endl;
	}
	cout <<     "*********************************************" << endl;
	cout << endl;
    }
       
    cout << "Done testing tstSphyramid_Mesh_Pack.cc on " << C4::node() << endl;
    
    C4::Finalize();
}   

//---------------------------------------------------------------------------//
//                        end of tstSphyramid_Mesh_Pack.cc.cc
//---------------------------------------------------------------------------//
