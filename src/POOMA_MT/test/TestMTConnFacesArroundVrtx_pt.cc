//----------------------------------*-C++-*----------------------------------//
// TestMTConnFacesArroundVrtx_pt.cc
// Randy M. Roberts
// Sat Aug 21 15:49:22 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// instantiate the template

#include "PoomaMesh_XYZFactory.hh"
#include "meshTest/TestMTConnFacesArroundVrtx.t.hh"

using rtt_meshTest::TestMTConnFacesArroundVrtx;
using rtt_POOMA_MT_test::PoomaMesh_XYZFactory;

typedef TestMTConnFacesArroundVrtx<PoomaMesh_XYZFactory> Tester;

// Only instantiate the run method.
// If you instantiate the entire class, methods that are not required
// for this particular mesh will be instantiated, causing unnecessary
// compiler errors.

template void Tester::run();

//---------------------------------------------------------------------------//
//                              end of TestMTConnFacesArroundVrtx_pt.cc
//---------------------------------------------------------------------------//
