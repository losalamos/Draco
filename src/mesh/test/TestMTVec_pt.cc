//----------------------------------*-C++-*----------------------------------//
// TestMTVec_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "Mesh_XYZFactory.hh"
#include "meshTest/TestMTVec.t.hh"

using rtt_meshTest::TestMTVec;
using rtt_mesh_test::Mesh_XYZFactory;

template class TestMTVec<Mesh_XYZFactory>;

//---------------------------------------------------------------------------//
//                              end of TestMTVec_pt.cc
//---------------------------------------------------------------------------//
