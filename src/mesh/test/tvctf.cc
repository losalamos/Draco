//----------------------------------*-C++-*----------------------------------//
// tvctf.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the Mesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "meshTest/TestMTFields.hh"
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
using rtt_meshTest::TestMTFields;

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
	    cout << "Initiating test of the MT::vctf."
		 << endl;

	    typedef TestMTFields<Mesh_XYZFactory> Tester;
	    
	    Tester tester(mtFactory, cout);
	    tester.run<Tester::VCTF>("vctf");

	    passed = tester.passed();
	}
    }
    catch( dsxx::assertion& a )
    {
	cout << "Failed assertion: " << a.what() << endl;
    }

    // Print the status of the test.

    cout << endl;
    cout <<     "******************************************" << endl;
    if (passed) 
    {
        cout << "**** Mesh_XYZ::vctf Self Test: PASSED ****" << endl;
    }
    else
    {
        cout << "**** Mesh_XYZ::vctf Self Test: FAILED ****" << endl;
    }
    cout <<     "******************************************" << endl;
    cout << endl;

    cout << "Done testing Mesh_XYZ::vctf.\n";

    C4::Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tvctf.cc
//---------------------------------------------------------------------------//
