//----------------------------------*-C++-*----------------------------------//
// tvctf.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the PoomaMesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "../PoomaMesh_XYZ.hh"
#include "PoomaMesh_XYZFactory.hh"

#include "meshTest/TestMTFields.hh"

#include "Tester.hh"

using rtt_meshTest::TestMTFields;

int main( int argc, char *argv[] )
{
    using namespace rtt_POOMA_MT_test;
    
    typedef PoomaMesh_XYZFactory Factory;
    typedef TestMTFields<Factory> Test;
    
    try
    {
	Tester<Factory> tester("test.in", argc, argv);
	tester.run<Test,Test::VCTF>("PoomaMesh_XYZ::vctf");
    }
    catch( dsxx::assertion& a )
    {
	std::cout << "Failed assertion: " << a.what() << std::endl;
    }
    
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tvctf.cc
//---------------------------------------------------------------------------//
