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

#include "IMC_Test.hh"
#include "IMC_DD_Test.hh"
#include "mc/Topology.hh"
#include "mc/OS_Mesh.hh"
#include "ds++/SP.hh"
#include "../Particle.hh"
#include "../Source_Builder.t.hh"
#include "../Rep_Source_Builder.t.hh"
#include "../DD_Source_Builder.t.hh"
#include "../Source.t.hh"
#include "../Source_Init.t.hh"

namespace rtt_imc
{

typedef rtt_dsxx::SP<rtt_imc_test::IMC_Flat_Interface>  SP_REP_IT;
typedef rtt_dsxx::SP<rtt_imc_dd_test::IMC_DD_Interface> SP_DD_IT;
typedef rtt_mc::OS_Mesh                                 MT;
typedef rtt_dsxx::SP<MT>                                SP_MT;
typedef rtt_dsxx::SP<rtt_mc::Topology>                  SP_TOP;
typedef Particle<MT>                                    PT;

template class Source_Builder<MT,PT>;
template Source_Builder<MT,PT>::Source_Builder(SP_REP_IT, SP_MT, SP_TOP);

template class Rep_Source_Builder<MT,PT>;
template Rep_Source_Builder<MT,PT>::Rep_Source_Builder(SP_REP_IT, SP_MT,
						       SP_TOP);

template class DD_Source_Builder<MT,PT>;
template DD_Source_Builder<MT,PT>::DD_Source_Builder(SP_DD_IT, SP_MT,
						     SP_TOP);

template class Source<MT,PT>;

template class Source_Init<rtt_imc_test::IMC_Flat_Interface, MT>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Source_pt.cc
//---------------------------------------------------------------------------//
