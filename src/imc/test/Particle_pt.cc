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
#include "rng/Rnd_Control.hh"
#include "../Particle.t.hh"

namespace rtt_imc
{

typedef rtt_mc::OS_Mesh       MT;
typedef rtt_imc::Particle<MT> PT;
typedef rtt_rng::Rnd_Control  RC;

template class Particle<MT>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Particle_pt.cc
//---------------------------------------------------------------------------//
