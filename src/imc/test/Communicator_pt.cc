//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/Communicator_pt.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 18 16:18:36 2000
 * \brief  Communicator class instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "../Particle.hh"
#include "../Communicator.t.hh"
#include "../Comm_Builder.t.hh"

namespace rtt_imc
{

typedef rtt_imc::Particle<rtt_mc::OS_Mesh> PT;

template class Communicator<PT>;

template class Comm_Builder<PT>;

} // end namespace rtt_imc


//---------------------------------------------------------------------------//
//                              end of Communicator_pt.cc
//---------------------------------------------------------------------------//
