//----------------------------------*-C++-*----------------------------------//
// TestMTFields_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "Mesh_XYZFactory.hh"
#include "meshTest/TestMTFields.t.hh"

using rtt_meshTest::TestMTFields;
using rtt_mesh_test::Mesh_XYZFactory;

template class TestMTFields<Mesh_XYZFactory>;

typedef TestMTFields<Mesh_XYZFactory> Tester;

template void Tester::run<Tester::BSTF>(const std::string&);
template void Tester::run<Tester::CCTF>(const std::string&);
template void Tester::run<Tester::FCDTF>(const std::string&);
template void Tester::run<Tester::NCTF>(const std::string&);
template void Tester::run<Tester::VCTF>(const std::string&);

//---------------------------------------------------------------------------//
//                              end of TestMTFields_pt.cc
//---------------------------------------------------------------------------//
