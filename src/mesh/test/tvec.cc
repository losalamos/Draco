//----------------------------------*-C++-*----------------------------------//
// tvec.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the Mesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "meshTest/TestMTVec.hh"

#include "../Mesh_XYZ.hh"
#include "Mesh_XYZFactory.hh"

#include "Tester.hh"

using rtt_meshTest::TestMTVec;

int main( int argc, char *argv[] )
{
    using namespace rtt_mesh_test;
    
    typedef Mesh_XYZFactory Factory;
    typedef TestMTVec<Factory> Test;
    
    try
    {
	Tester<Factory> tester("test.in", argc, argv);
	tester.run<Test>("Mesh_XYZ::ccvsf::value_type container");
    }
    catch( rtt_dsxx::assertion& a )
    {
	std::cout << "Failed assertion: " << a.what() << std::endl;
    }
    
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tvec.cc
//---------------------------------------------------------------------------//
