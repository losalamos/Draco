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
#include "mc/RZWedge_Mesh.hh"
#include "../Gray_Particle.t.hh"
#include "../Multigroup_Particle.t.hh"

namespace rtt_imc
{

typedef rtt_mc::OS_Mesh      OS_MT;
typedef rtt_mc::RZWedge_Mesh RZ_MT;

template class Gray_Particle<OS_MT>;
template class Multigroup_Particle<OS_MT>;

template class Gray_Particle<RZ_MT>;
template class Multigroup_Particle<RZ_MT>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Particle_pt.cc
//---------------------------------------------------------------------------//
