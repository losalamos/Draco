//----------------------------------*-C++-*----------------------------------//
// tcctf.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the Mesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "meshTest/TestMTFields.hh"
#include "../Mesh_XYZ.hh"
#include "Mesh_XYZFactory.hh"

#include "Tester.hh"

using rtt_meshTest::TestMTFields;

int main( int argc, char *argv[] )
{
    using namespace rtt_mesh_test;
    
    typedef Mesh_XYZFactory Factory;
    typedef TestMTFields<Factory> Test;
    
    try
    {
	Tester<Factory> tester("test.in", argc, argv);
	tester.run<Test,Test::CCTF>("Mesh_XYZ::cctf");
    }
    catch( dsxx::assertion& a )
    {
	std::cout << "Failed assertion: " << a.what() << std::endl;
    }
    
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tcctf.cc
//---------------------------------------------------------------------------//
