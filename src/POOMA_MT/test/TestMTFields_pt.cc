//----------------------------------*-C++-*----------------------------------//
// TestMTFields_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "PoomaMesh_XYZFactory.hh"
#include "meshTest/TestMTFields.t.hh"

using rtt_meshTest::TestMTFields;
using rtt_POOMA_MT_test::PoomaMesh_XYZFactory;

template class TestMTFields<PoomaMesh_XYZFactory>;

typedef TestMTFields<PoomaMesh_XYZFactory> Tester;

template void Tester::run<Tester::BSTF>();
template void Tester::run<Tester::CCTF>();
template void Tester::run<Tester::FCDTF>();
template void Tester::run<Tester::NCTF>();
template void Tester::run<Tester::VCTF>();

//---------------------------------------------------------------------------//
//                              end of TestMTFields_pt.cc
//---------------------------------------------------------------------------//
