//----------------------------------*-C++-*----------------------------------//
// Communicator_pt.cc
// Thomas M. Evans
// Mon Jul 13 10:35:24 1998
//---------------------------------------------------------------------------//
// @> Instantialyzer for template class Communicator
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "Particle.hh"
#include "Communicator.t.hh"

namespace rtt_imc 
{

template class Communicator<Particle<OS_Mesh> >;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Communicator_pt.cc
//---------------------------------------------------------------------------//
