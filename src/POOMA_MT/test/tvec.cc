//----------------------------------*-C++-*----------------------------------//
// tvec.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the PoomaMesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "../PoomaMesh_XYZ.hh"
#include "PoomaMesh_XYZFactory.hh"

#include "meshTest/TestMTVec.hh"

#include "Tester.hh"

using rtt_meshTest::TestMTVec;

int main( int argc, char *argv[] )
{
    using namespace rtt_POOMA_MT_test;
    
    typedef PoomaMesh_XYZFactory Factory;
    typedef TestMTVec<Factory> Test;
    
    try
    {
	Tester<Factory> tester("test.in", argc, argv);
	tester.run<Test>("PoomaMesh_XYZ::ccvsf::value_type container");
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
