//----------------------------------*-C++-*----------------------------------//
// FieldTester_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "Mesh_XYZFactory.hh"
#include "meshTest/TestMTFields.hh"
#include "meshTest/FieldTester.t.hh"

using rtt_meshTest::TestMTFields;
using rtt_meshTest::FieldTester;
using rtt_mesh_test::Mesh_XYZFactory;

typedef TestMTFields<Mesh_XYZFactory> Tester;

template void FieldTester<Mesh_XYZFactory::MT, Tester::BSTF>::run();
template void FieldTester<Mesh_XYZFactory::MT, Tester::CCTF>::run();
template void FieldTester<Mesh_XYZFactory::MT, Tester::FCDTF>::run();
template void FieldTester<Mesh_XYZFactory::MT, Tester::NCTF>::run();
template void FieldTester<Mesh_XYZFactory::MT, Tester::VCTF>::run();

//---------------------------------------------------------------------------//
//                              end of FieldTester_pt.cc
//---------------------------------------------------------------------------//
