//----------------------------------*-C++-*----------------------------------//
// TestMTModel_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "Mesh_XYZFactory.hh"
#include "meshTest/TestMTModel.t.hh"

using rtt_meshTest::TestMTModel;
using rtt_mesh_test::Mesh_XYZFactory;

template void TestMTModel<Mesh_XYZFactory>::run();

//---------------------------------------------------------------------------//
//                              end of TestMTModel_pt.cc
//---------------------------------------------------------------------------//
