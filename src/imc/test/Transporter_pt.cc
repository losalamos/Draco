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
#include "../Gray_Particle.hh"
#include "../Multigroup_Particle.hh"
#include "../Frequency.hh"
#include "../Rep_Transporter.t.hh"
#include "../DD_Transporter.t.hh"

namespace rtt_imc
{

typedef rtt_mc::OS_Mesh         MT;
typedef Gray_Frequency          Gray;
typedef Multigroup_Frequency    MG;
typedef Gray_Particle<MT>       GPT;
typedef Multigroup_Particle<MT> MGPT;

template class Rep_Transporter<MT,Gray,GPT>;
template class DD_Transporter<MT,Gray,GPT>;

template class Rep_Transporter<MT,MG,MGPT>;
template class DD_Transporter<MT,MG,MGPT>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Transporter_pt.cc
//---------------------------------------------------------------------------//
