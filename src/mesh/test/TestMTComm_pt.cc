//----------------------------------*-C++-*----------------------------------//
// TestMTComm_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "Mesh_XYZFactory.hh"
#include "meshTest/TestMTComm.t.hh"

using rtt_meshTest::TestMTComm;
using rtt_mesh_test::Mesh_XYZFactory;

template class TestMTComm<Mesh_XYZFactory>;

//---------------------------------------------------------------------------//
//                              end of TestMTComm_pt.cc
//---------------------------------------------------------------------------//
