//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstSource_Builder_pt.cc
 * \author Thomas M. Evans
 * \date   Wed Dec  8 16:43:40 1999
 * \brief  Source_Builder explicit instantiations for testing.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "IMC_Test.hh"
#include "IMC_DD_Test.hh"
#include "mc/Topology.hh"
#include "mc/OS_Mesh.hh"
#include "ds++/SP.hh"
#include "../Source_Builder.t.hh"
#include "../Rep_Source_Builder.t.hh"
#include "../DD_Source_Builder.t.hh"
#include "../Source.t.hh"
#include "../Source_Init.t.hh"

namespace rtt_imc
{

using rtt_imc_test::IMC_Interface;
using rtt_imc_dd_test::IMC_DD_Interface;
using rtt_mc::OS_Mesh;
using rtt_mc::Topology;
using rtt_dsxx::SP;

template class Source_Builder<OS_Mesh>;
template Source_Builder<OS_Mesh>::Source_Builder(SP<IMC_Interface>, 
						 SP<OS_Mesh>, 
						 SP<Topology>);

template class Rep_Source_Builder<OS_Mesh>;
template Rep_Source_Builder<OS_Mesh>::Rep_Source_Builder(SP<IMC_Interface>,
							 SP<OS_Mesh>,
							 SP<Topology>);

template class DD_Source_Builder<OS_Mesh>;
template DD_Source_Builder<OS_Mesh>::DD_Source_Builder(SP<IMC_DD_Interface>,
						       SP<OS_Mesh>,
						       SP<Topology>);

template class Source<OS_Mesh>;

template class Source_Init<IMC_Interface, OS_Mesh>;

} // end of namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of tstSource_Builder_pt.cc
//---------------------------------------------------------------------------//
