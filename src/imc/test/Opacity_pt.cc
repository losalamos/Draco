//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/Opacity_pt.cc
 * \author Thomas M. Evans
 * \date   Fri Aug 18 16:21:43 2000
 * \brief  Opacity class instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "mc/RZWedge_Mesh.hh"
#include "../Opacity_Builder.t.hh"
#include "../Opacity.t.hh"
#include "../Mat_State.t.hh"
#include <iostream> 

namespace rtt_imc
{

typedef rtt_mc::OS_Mesh      OS_MT;
typedef rtt_mc::RZWedge_Mesh RZ_MT;

using std::ostream;

// OS_Mesh instantiations

template class Opacity_Builder<OS_MT>;

template class Opacity<OS_MT>;
template ostream& operator<<(ostream &, const Opacity<OS_MT> &);

template class Mat_State<OS_MT>;
template ostream& operator<<(ostream &, const Mat_State<OS_MT> &);

// RZWedge_Mesh instantiations

template class Mat_State<RZ_MT>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Opacity_pt.cc
//---------------------------------------------------------------------------//
