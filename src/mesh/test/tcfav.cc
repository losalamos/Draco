//----------------------------------*-C++-*----------------------------------//
// tcav.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the Mesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "meshTest/TestMTConnFacesArroundVrtx.hh"
#include "../Mesh_XYZ.hh"
#include "Mesh_XYZFactory.hh"

#include "Tester.hh"

using rtt_meshTest::TestMTConnFacesArroundVrtx;

int main( int argc, char *argv[] )
{
    using namespace rtt_mesh_test;
    
    typedef Mesh_XYZFactory Factory;
    typedef TestMTConnFacesArroundVrtx<Factory> Test;
    
    try
    {
	Tester<Factory> tester("test.in", argc, argv);
	tester.run<Test>("Mesh_XYZ connection faces arround vertices");
    }
    catch( dsxx::assertion& a )
    {
	std::cout << "Failed assertion: " << a.what() << std::endl;
    }
    
    return 0;
}

//---------------------------------------------------------------------------//
//                              end of tcav.cc
//---------------------------------------------------------------------------//
