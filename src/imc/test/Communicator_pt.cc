//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/Communicator_pt.cc
 * \author Thomas M. Evans
 * \date   Tue Jan  8 13:53:42 2002
 * \brief  Communicator and Particle_Buffer instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../Gray_Particle.hh"
#include "../Multigroup_Particle.hh"
#include "mc/OS_Mesh.hh"
#include "mc/Communicator.t.hh"
#include "mc/Particle_Buffer.t.hh"

namespace rtt_mc
{

typedef rtt_imc::Gray_Particle<OS_Mesh> PT;

template class Communicator<PT>;
template class Particle_Buffer<PT>;
template class Recv_Particle_Buffer<PT>;
template class Send_Particle_Buffer<PT>;

typedef rtt_imc::Multigroup_Particle<OS_Mesh> MGPT;

template class Communicator<MGPT>;
template class Particle_Buffer<MGPT>;
template class Recv_Particle_Buffer<MGPT>;
template class Send_Particle_Buffer<MGPT>;

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of Communicator_pt.cc
//---------------------------------------------------------------------------//
