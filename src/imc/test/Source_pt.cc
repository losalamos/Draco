//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/Source_pt.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 18 16:26:04 2000
 * \brief  Source class instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "mc/RZWedge_Mesh.hh"
#include "../Gray_Particle.hh"
#include "../Multigroup_Particle.hh"
#include "../Frequency.hh"
#include "../Source_Builder.t.hh"
#include "../Rep_Source_Builder.t.hh"
#include "../DD_Source_Builder.t.hh"
#include "../Source.t.hh"

namespace rtt_imc
{

typedef rtt_mc::OS_Mesh         MT;
typedef rtt_mc::RZWedge_Mesh    RZ;
typedef Gray_Frequency          G;
typedef Multigroup_Frequency    MG;
typedef Gray_Particle<MT>       GPT;
typedef Multigroup_Particle<MT> MGPT;
typedef Multigroup_Particle<RZ> MGPTRZ;

template class Source_Builder<MT,G,GPT>;
template class Rep_Source_Builder<MT,G,GPT>;
template class DD_Source_Builder<MT,G,GPT>;

template class Source_Builder<MT,MG,MGPT>;
template class Rep_Source_Builder<MT,MG,MGPT>;
template class DD_Source_Builder<MT,MG,MGPT>;

template class Source<MT,G,GPT>;
template class Source<MT,MG,MGPT>;

template class Source_Builder<RZ,MG,MGPTRZ>;
template class Rep_Source_Builder<RZ,MG,MGPTRZ>;
template class DD_Source_Builder<RZ,MG,MGPTRZ>;
template class Source<RZ,MG,MGPTRZ>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Source_pt.cc
//---------------------------------------------------------------------------//
