//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   imc/test/Particle_pt.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 18 16:14:04 2000
 * \brief  Particle class instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "../Gray_Particle.t.hh"
#include "../Multigroup_Particle.t.hh"

namespace rtt_imc
{

typedef rtt_mc::OS_Mesh MT;

template class Gray_Particle<MT>;

template class Multigroup_Particle<MT>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Particle_pt.cc
//---------------------------------------------------------------------------//
