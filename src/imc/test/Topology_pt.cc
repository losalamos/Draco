//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/Topology_pt.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 18 17:35:48 2000
 * \brief  Topology class instantiations.
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

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Topology_pt.cc
//---------------------------------------------------------------------------//
