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

#include "../Mesh_XYZ.hh"
#include "Mesh_XYZFactory.hh"

#include "Tester.hh"

using rtt_meshTest::TestMTModel;
using rtt_meshTest::TestMTComm;

int main( int argc, char *argv[] )
{
    using namespace rtt_mesh_test;
    
    typedef Mesh_XYZFactory Factory;
    typedef TestMTModel<Factory> Test1;
    typedef TestMTComm<Factory> Test2;
    
    Tester<Factory> tester("test.in", argc, argv);
    tester.run<Test1>("Mesh_XYZ Model of MT");
    tester.run<Test2>("Mesh_XYZ communication");

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of txyz.cc
//---------------------------------------------------------------------------//
