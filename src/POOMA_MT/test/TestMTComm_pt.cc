//----------------------------------*-C++-*----------------------------------//
// TestMTComm_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "PoomaMesh_XYZFactory.hh"
#include "meshTest/TestMTComm.t.hh"

using rtt_meshTest::TestMTComm;
using rtt_POOMA_MT_test::PoomaMesh_XYZFactory;

template void TestMTComm<PoomaMesh_XYZFactory>::run();

//---------------------------------------------------------------------------//
//                              end of TestMTComm_pt.cc
//---------------------------------------------------------------------------//
