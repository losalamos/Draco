//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/Transporter_pt.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 18 17:37:20 2000
 * \brief  Transporter class instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "../Particle.hh"
#include "../Rep_Transporter.t.hh"
#include "../DD_Transporter.t.hh"

namespace rtt_imc
{

typedef rtt_mc::OS_Mesh MT;
typedef Particle<MT>    PT;

template class Rep_Transporter<MT, PT>;

template class DD_Transporter<MT, PT>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Transporter_pt.cc
//---------------------------------------------------------------------------//
