//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/DD_Mesh.hh
 * \author Thomas M. Evans
 * \date   Thu May  4 14:57:15 2000
 * \brief  Build a DD Mesh for testing.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __test_DD_Mesh_hh__
#define __test_DD_Mesh_hh__

#include "../OS_Mesh.hh"
#include "../General_Topology.hh"
#include "ds++/SP.hh"

namespace rtt_mc_test
{
 
//===========================================================================//
// BUILD DD MESHES
//===========================================================================//
// build 9 cell meshes 2 cells on proc 0; 2 cells on proc 1; 2 cells on proc
// 2; 3 cells on proc 3

rtt_dsxx::SP<rtt_mc::OS_Mesh> build_Mesh();


//===========================================================================//
// BUILD DD TOPOLOGIES
//===========================================================================//
// build a General Topology for 9 cell mesh described above

rtt_dsxx::SP<rtt_mc::Topology> build_Topology();

} // end namespace rtt_mc_test

#endif                          // __test_DD_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of test/DD_Mesh.hh
//---------------------------------------------------------------------------//
