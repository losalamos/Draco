//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstTransporter_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Apr 17 16:13:12 2000
 * \brief  Transporter instantiations for testing.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "../Particle.hh"
#include "../Rep_Transporter.t.hh"
#include "../DD_Transporter.t.hh"
#include "../Communicator.t.hh"
#include "../Source.t.hh"

namespace rtt_imc
{

using rtt_mc::OS_Mesh;

template class Rep_Transporter<OS_Mesh, Particle<OS_Mesh> >;

template class DD_Transporter<OS_Mesh, Particle<OS_Mesh> >;

template class Communicator<Particle<OS_Mesh> >;

template class Source<OS_Mesh, Particle<OS_Mesh> >;

} // end of rtt_imc

//---------------------------------------------------------------------------//
//                              end of tstTransporter_pt.cc
//---------------------------------------------------------------------------//
