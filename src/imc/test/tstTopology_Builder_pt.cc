//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstTopology_Builder_pt.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 30 10:37:32 1999
 * \brief  Topology_Builder explicit instantiations for testing.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "../Topology_Builder.t.hh"

namespace rtt_imc
{

using rtt_mc::OS_Mesh;

template class Topology_Builder<OS_Mesh>;

} // end of namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of tstTopology_Builder_pt.cc
//---------------------------------------------------------------------------//
