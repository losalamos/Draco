//----------------------------------*-C++-*----------------------------------//
// Particle_Buffer_pt.cc
// Thomas M. Evans
// Tue May 12 16:00:54 1998
//---------------------------------------------------------------------------//
// @> Instantialyzer for template class Particle_Buffer
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "Particle.hh"
#include "Particle_Buffer.t.hh"

IMCSPACE

using rtt_rng::Rnd_Control;

template class Particle_Buffer<Particle<OS_Mesh> >;

template
Particle_Buffer<Particle<OS_Mesh> >::Particle_Buffer(const OS_Mesh &, 
						     const Rnd_Control &);  

CSPACE

//---------------------------------------------------------------------------//
//                              end of Particle_Buffer_pt.cc
//---------------------------------------------------------------------------//
