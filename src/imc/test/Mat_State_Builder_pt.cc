//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/Mat_State_Builder_pt.cc
 * \author Thomas M. Evans
 * \date   Tue Jan 22 14:54:02 2002
 * \brief  Mat_State_Builder instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "../Frequency.hh"
#include "../CDI_Mat_State_Builder.t.hh"
#include "../Flat_Mat_State_Builder.t.hh"

namespace rtt_imc
{

typedef rtt_mc::OS_Mesh MT;

template class CDI_Mat_State_Builder<MT, Gray_Frequency>;

template class CDI_Mat_State_Builder<MT, Multigroup_Frequency>;

template class Flat_Mat_State_Builder<MT, Gray_Frequency>;

template class Flat_Mat_State_Builder<MT, Multigroup_Frequency>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Mat_State_Builder_pt.cc
//---------------------------------------------------------------------------//
