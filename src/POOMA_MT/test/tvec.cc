//----------------------------------*-C++-*----------------------------------//
// tvec.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the PoomaMesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "meshTest/TestMTVec.hh"

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
using rtt_meshTest::TestMTVec;

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
	    cout << "Initiating test of the "
		 << "PoomaMesh_XYZ::ccvsf::value_type container."
		 << endl;

	    TestMTVec<PoomaMesh_XYZFactory> tester(mtFactory, cout);
	    tester.run();

	    passed = tester.passed();
	}

    }
    catch( dsxx::assertion& a )
    {
	cout << "Failed assertion: " << a.what() << endl;
    }

    // Print the status of the test.

    cout << endl;
    cout << "***********************************************************"
         << endl;
    if (passed) 
    {
        cout <<
            "**** MT::ccvsf::value_type Container Self Test: PASSED ****"
             << endl;
    }
    else
    {
        cout <<
            "**** MT::ccvsf::value_type Container Self Test: FAILED ****"
             << endl;
    }
    cout << "***********************************************************"
         << endl;
    cout << endl;

    cout << "Done testing MT::ccvsf::value_type container.\n";

    C4::Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tvec.cc
//---------------------------------------------------------------------------//
