//----------------------------------*-C++-*----------------------------------//
// tstParticle_pt.cc
// Thomas M. Evans
// Thu Apr 29 15:35:23 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Particle and Particle_Buffer instantiations
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "rng/Rnd_Control.hh"
#include "../Particle.t.hh"
#include "../Particle_Buffer.t.hh"

namespace rtt_imc
{

using rtt_mc::OS_Mesh;
using rtt_rng::Rnd_Control;

// particle class instantiation

template class Particle<OS_Mesh>;

// random number class instantiation

template class Particle_Buffer<Particle<OS_Mesh> >;

template
Particle_Buffer<Particle<OS_Mesh> >::Particle_Buffer(const OS_Mesh &, 
						     const Rnd_Control &);
}  // end of namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of tstParticle_pt.cc
//---------------------------------------------------------------------------//
