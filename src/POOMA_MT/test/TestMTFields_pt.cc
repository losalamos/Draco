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

template void Tester::run<Tester::BSTF>(const std::string&);
template void Tester::run<Tester::CCTF>(const std::string&);
template void Tester::run<Tester::FCDTF>(const std::string&);
template void Tester::run<Tester::NCTF>(const std::string&);
template void Tester::run<Tester::VCTF>(const std::string&);

//---------------------------------------------------------------------------//
//                              end of TestMTFields_pt.cc
//---------------------------------------------------------------------------//
