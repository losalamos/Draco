//----------------------------------*-C++-*----------------------------------//
// txyz.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the Mesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "meshTest/TestMTModel.hh"
#include "meshTest/TestMTComm.hh"

#include "meshTest/Release.hh"
#include "../Release.hh"
#include "TestRelease.hh"

#include "../Mesh_XYZ.hh"
#include "Mesh_XYZFactory.hh"

#include "nml/Group.hh"

#include "c4/global.hh"

#include <iostream>
using std::cout;
using std::endl;

using rtt_mesh_test::Mesh_XYZFactory;
using rtt_meshTest::TestMTModel;
using rtt_meshTest::TestMTComm;

void version(const std::string &progname)
{
    cout << progname << ": version " << rtt_mesh_test::release() << endl;
    cout << "mesh: version " << rtt_mesh::release() << endl;
    cout << "meshTest: version " << rtt_meshTest::release() << endl;
}

int main( int argc, char *argv[] )
{
    C4::Init( argc, argv );

    for (int arg=1; arg < argc; arg++)
    {
	if (std::string(argv[arg]) == "--version")
	{
	    version(argv[0]);
	    C4::Finalize();
	    return 0;
	}
    }
    
    bool passed = false;
    
    try {

	NML_Group g( "test" );
	Mesh_DB mdb;
	mdb.setup_namelist( g );
	g.readgroup( "test.in" );

	Mesh_XYZFactory mtFactory(mdb);

	{
	    cout << "Initiating test of the Mesh_XYZ model of the MT concept."
		 << endl;

	    TestMTModel<Mesh_XYZFactory> tester(mtFactory, cout);
	    tester.run();

	    passed = tester.passed();
	}

	{
	
	    cout << "Initiating test of the Mesh_XYZ communication."
		 << endl;

	    TestMTComm<Mesh_XYZFactory> tester(mtFactory, cout);
	    tester.run();

	    passed = passed && tester.passed();
	}
	
    }
    catch( dsxx::assertion& a )
    {
	cout << "Failed assertion: " << a.what() << endl;
    }

    // Print the status of the test.

    cout << endl;
    cout <<     "************************************" << endl;
    if (passed) 
    {
        cout << "**** Mesh_XYZ Self Test: PASSED ****" << endl;
    }
    else
    {
        cout << "**** Mesh_XYZ Self Test: FAILED ****" << endl;
    }
    cout <<     "************************************" << endl;
    cout << endl;

    cout << "Done testing Mesh_XYZ.\n";

    C4::Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of txyz.cc
//---------------------------------------------------------------------------//
