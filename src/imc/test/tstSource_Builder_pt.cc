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

#include "mc/OS_Mesh.hh"
#include "../Source_Builder.t.hh"

namespace rtt_imc
{

using rtt_mc::OS_Mesh;

template class Source_Builder<OS_Mesh>;

} // end of namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of tstSource_Builder_pt.cc
//---------------------------------------------------------------------------//
