//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/Random_Walk_pt.cc
 * \author Thomas M. Evans
 * \date   Tue Jan 14 18:50:01 2003
 * \brief  Random_Walk instantiation.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../Random_Walk.t.hh"
#include "mc/OS_Mesh.hh"
#include "mc/RZWedge_Mesh.hh"

namespace rtt_imc
{

typedef rtt_mc::OS_Mesh      MT;
typedef rtt_mc::RZWedge_Mesh RZMT;

template class Random_Walk<MT>;
template class Random_Walk<RZMT>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Random_Walk_pt.cc
//---------------------------------------------------------------------------//
