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


//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void test_pack()
{
    using rtt_mc_test::Parser;
    using rtt_dsxx::SP;
    using rtt_mc::Sphyramid_Mesh;
    using rtt_mc::Sphyramid_Builder;

    // make a builder from parsing the Sphyramid_Input file
    SP<Parser> parser(new Parser("Sphyramid_Input"));
    Sphyramid_Builder builder(parser);

    SP<Sphyramid_Mesh> m1;
    SP<Sphyramid_Mesh> m2;
    SP<Sphyramid_Mesh> m3;

    // build the mesh
    m1 = builder.build_Mesh();
    if (m1->num_cells() != 5) ITFAILS;

    // replicate the mesh
    {
	SP<Sphyramid_Mesh::Pack> pack = m1->pack();
	if (pack->get_num_packed_cells() != 5) ITFAILS;
	//m2 = pack->unpack();
    }


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
