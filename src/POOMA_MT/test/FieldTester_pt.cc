//----------------------------------*-C++-*----------------------------------//
// FieldTester_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "PoomaMesh_XYZFactory.hh"
#include "meshTest/TestMTFields.hh"
#include "meshTest/FieldTester.t.hh"

using rtt_meshTest::TestMTFields;
using rtt_meshTest::FieldTester;
using rtt_POOMA_MT_test::PoomaMesh_XYZFactory;

typedef TestMTFields<PoomaMesh_XYZFactory> Tester;

template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::BSTF>::run();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::CCTF>::run();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::FCDTF>::run();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::NCTF>::run();
template void FieldTester<PoomaMesh_XYZFactory::MT, Tester::VCTF>::run();

//---------------------------------------------------------------------------//
//                              end of FieldTester_pt.cc
//---------------------------------------------------------------------------//
