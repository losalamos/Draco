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

#include "IMC_Test.hh"
#include "mc/OS_Mesh.hh"
#include "ds++/SP.hh"
#include "../Topology_Builder.t.hh"

namespace rtt_imc
{

using rtt_imc_test::IMC_Interface;
using rtt_mc::OS_Mesh;
using dsxx::SP;

template class Topology_Builder<OS_Mesh>;

template Topology_Builder<OS_Mesh>::Topology_Builder(SP<IMC_Interface>);

} // end of namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of tstTopology_Builder_pt.cc
//---------------------------------------------------------------------------//
