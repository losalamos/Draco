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
#include "mc/Topology.hh"
#include "mc/OS_Mesh.hh"
#include "ds++/SP.hh"
#include "../Source_Builder.t.hh"
#include "../Rep_Source_Builder.t.hh"
#include "../Source.t.hh"

namespace rtt_imc
{

using rtt_imc_test::IMC_Interface;
using rtt_mc::OS_Mesh;
using rtt_mc::Topology;
using dsxx::SP;

template class Source_Builder<OS_Mesh>;
template Source_Builder<OS_Mesh>::Source_Builder(SP<IMC_Interface>, 
						 SP<OS_Mesh>, 
						 SP<Topology>);

template class Rep_Source_Builder<OS_Mesh>;
template Rep_Source_Builder<OS_Mesh>::Rep_Source_Builder(SP<IMC_Interface>,
							 SP<OS_Mesh>,
							 SP<Topology>);

template class Source<OS_Mesh>;

} // end of namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of tstSource_Builder_pt.cc
//---------------------------------------------------------------------------//
