//----------------------------------*-C++-*----------------------------------//
// tvctf.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the PoomaMesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "meshTest/TestMTFields.hh"

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
using rtt_meshTest::TestMTFields;

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
	    cout << "Initiating test of the MT::vctf."
		 << endl;

	    typedef TestMTFields<PoomaMesh_XYZFactory> Tester;
	    
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
        cout << "**** PoomaMesh_XYZ::vctf Self Test: PASSED ****" << endl;
    }
    else
    {
        cout << "**** PoomaMesh_XYZ::vctf Self Test: FAILED ****" << endl;
    }
    cout <<     "******************************************" << endl;
    cout << endl;

    cout << "Done testing PoomaMesh_XYZ::vctf.\n";

    C4::Finalize();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tvctf.cc
//---------------------------------------------------------------------------//
