//----------------------------------*-C++-*----------------------------------//
// txyz.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the PoomaMesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "meshTest/TestMTModel.hh"
#include "meshTest/TestMTComm.hh"

#include "../PoomaMesh_XYZ.hh"
#include "PoomaMesh_XYZFactory.hh"

#include "nml/Group.hh"

#include "c4/global.hh"

#include <iostream>
using std::cout;
using std::endl;

#include "utils.hh"
using rtt_POOMA_MT_test::getMTFactory;
using rtt_POOMA_MT_test::version;

using rtt_POOMA_MT_test::PoomaMesh_XYZFactory;
using rtt_meshTest::TestMTModel;
using rtt_meshTest::TestMTComm;

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

	PoomaMesh_XYZFactory mtFactory = getMTFactory("test.in");

	{
	    cout << "Initiating test of the PoomaMesh_XYZ model"
		 << " of the MT concept."
		 << endl;

	    TestMTModel<PoomaMesh_XYZFactory> tester(mtFactory, cout);
	    tester.run();

	    passed = tester.passed();
	}

	{
	
	    cout << "Initiating test of the PoomaMesh_XYZ communication."
		 << endl;

	    TestMTComm<PoomaMesh_XYZFactory> tester(mtFactory, cout);
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
        cout << "**** PoomaMesh_XYZ Self Test: PASSED ****" << endl;
    }
    else
    {
        cout << "**** PoomaMesh_XYZ Self Test: FAILED ****" << endl;
    }
    cout <<     "************************************" << endl;
    cout << endl;

    cout << "Done testing PoomaMesh_XYZ.\n";

    C4::Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of txyz.cc
//---------------------------------------------------------------------------//
