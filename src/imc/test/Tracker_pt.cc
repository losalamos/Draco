//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/Tracker_pt.cc
 * \author Thomas M. Evans
 * \date   Wed Aug 20 15:54:18 2003
 * \brief  Explicit instantiation of surface trackers.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "mc/RZWedge_Mesh.hh"
#include "../Extrinsic_Tracker_Builder.t.hh"

namespace rtt_imc
{

template class Extrinsic_Tracker_Builder<rtt_mc::OS_Mesh>;

template class Extrinsic_Tracker_Builder<rtt_mc::RZWedge_Mesh>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                 end of Tracker_pt.cc
//---------------------------------------------------------------------------//
