//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstCommunicator_pt.cc
 * \author Thomas M. Evans
 * \date   Wed Jun  7 13:11:38 2000
 * \brief  tstCommunicator instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "rng/Rnd_Control.hh"
#include "../Particle.t.hh"
#include "../Comm_Builder.t.hh"
#include "../Communicator.t.hh"
#include "../Particle_Buffer.t.hh"
#include "../Tally.t.hh"

namespace rtt_imc
{

using rtt_mc::OS_Mesh;
using rtt_rng::Rnd_Control;

template class Communicator<Particle<OS_Mesh> >;

template class Comm_Builder<Particle<OS_Mesh> >;

// particle class instantiation

template class Particle<OS_Mesh>;

// particle buffer

template class Particle_Buffer<Particle<OS_Mesh> >;

template
Particle_Buffer<Particle<OS_Mesh> >::Particle_Buffer(const OS_Mesh &, 
						     const Rnd_Control &);

template class Tally<OS_Mesh>;

} // end namespace rtt_imc


//---------------------------------------------------------------------------//
//                              end of tstCommunicator_pt.cc
//---------------------------------------------------------------------------//
