//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/Opacity_pt.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 18 16:21:43 2000
 * \brief  Opacity class instantiations.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "mc/RZWedge_Mesh.hh"
#include "../Frequency.hh"
#include "../Opacity.t.hh"
#include "../Mat_State.t.hh"

namespace rtt_imc
{

typedef rtt_mc::OS_Mesh      OS_MT;
typedef rtt_mc::RZWedge_Mesh RZ_MT;

// OS_Mesh instantiations

template class Opacity<OS_MT, Gray_Frequency>;

template class Opacity<OS_MT, Multigroup_Frequency>; 

template class Mat_State<OS_MT>;

// RZWedge_Mesh instantiations

template class Opacity<RZ_MT, Gray_Frequency>;

template class Opacity<RZ_MT, Multigroup_Frequency>;

template class Mat_State<RZ_MT>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Opacity_pt.cc
//---------------------------------------------------------------------------//
