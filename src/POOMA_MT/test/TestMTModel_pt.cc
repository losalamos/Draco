//----------------------------------*-C++-*----------------------------------//
// TestMTModel_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "PoomaMesh_XYZFactory.hh"
#include "meshTest/TestMTModel.t.hh"

using rtt_meshTest::TestMTModel;
using rtt_POOMA_MT_test::PoomaMesh_XYZFactory;

template void TestMTModel<PoomaMesh_XYZFactory>::run();

//---------------------------------------------------------------------------//
//                              end of TestMTModel_pt.cc
//---------------------------------------------------------------------------//
