//----------------------------------*-C++-*----------------------------------//
// StructuredMeshCFAVTest_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "Mesh_XYZFactory.hh"
#include "meshTest/StructuredMeshCFAVTest.t.hh"

using rtt_meshTest::StructuredMeshCFAVTest;
using rtt_mesh_test::Mesh_XYZFactory;

typedef StructuredMeshCFAVTest<Mesh_XYZFactory> Tester;

// Only instantiate the run method.
// If you instantiate the entire class, methods that are not required
// for this particular mesh will be instantiated, causing unnecessary
// compiler errors.

template void Tester::run();

//---------------------------------------------------------------------------//
//                              end of StructuredMeshCFAVTest_pt.cc
//---------------------------------------------------------------------------//
