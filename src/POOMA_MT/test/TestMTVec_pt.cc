//----------------------------------*-C++-*----------------------------------//
// TestMTVec_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "PoomaMesh_XYZFactory.hh"
#include "meshTest/TestMTVec.t.hh"

using rtt_meshTest::TestMTVec;
using rtt_POOMA_MT_test::PoomaMesh_XYZFactory;

template void TestMTVec<PoomaMesh_XYZFactory>::run();

//---------------------------------------------------------------------------//
//                              end of TestMTVec_pt.cc
//---------------------------------------------------------------------------//
