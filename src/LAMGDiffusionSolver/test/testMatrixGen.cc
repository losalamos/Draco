//----------------------------------*-C++-*----------------------------------//
// testMatrixGen.cc
// Randy M. Roberts
// $Id$
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// This program tests the PoomaMesh_XYZ class as a model of MT.
//---------------------------------------------------------------------------//

#include "mesh/Mesh_XYZ.hh"
#include "Mesh_XYZFactory.hh"

#include "TestMatrixGen.hh"

int main( int argc, char *argv[] )
{
    using namespace rtt_LAMGDiffusionSolver_test;
    
    typedef Mesh_XYZFactory Factory;
    typedef TestMatrixGen<Factory> Tester;
    
    Tester tester("test.in", argc, argv, std::cout);
    tester.run();

    return 0;
}

//---------------------------------------------------------------------------//
//                              end of testMatrixGen.cc
//---------------------------------------------------------------------------//
