//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/test/PT_Instantiations_pt.cc
 * \author Thomas M. Evans
 * \date   Fri Dec 21 09:49:17 2001
 * \brief  Instantiations on PT.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "MC_Test.hh"
#include "../Particle_Buffer.t.hh"

namespace rtt_mc
{

typedef rtt_mc_test::Dummy_Particle PT;

template class Particle_Buffer<PT>;
template class Recv_Particle_Buffer<PT>;
template class Send_Particle_Buffer<PT>;

} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                              end of PT_Instantiations_pt.cc
//---------------------------------------------------------------------------//
