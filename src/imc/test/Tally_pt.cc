//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/Tally_pt.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 18 16:20:15 2000
 * \brief  Tally class instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "../Tally.t.hh"

namespace rtt_imc
{

using rtt_mc::OS_Mesh;

template class Tally<OS_Mesh>;

} // end namespace rtt_imc


//---------------------------------------------------------------------------//
//                              end of Tally_pt.cc
//---------------------------------------------------------------------------//
