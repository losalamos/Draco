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
#include "../Particle.hh"
#include "../Source_Builder.t.hh"
#include "../Rep_Source_Builder.t.hh"
#include "../DD_Source_Builder.t.hh"
#include "../Source.t.hh"

namespace rtt_imc
{

typedef rtt_mc::OS_Mesh MT;
typedef Particle<MT>    PT;

template class Source_Builder<MT,PT>;
template class Rep_Source_Builder<MT,PT>;
template class DD_Source_Builder<MT,PT>;

template class Source<MT,PT>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Source_pt.cc
//---------------------------------------------------------------------------//
