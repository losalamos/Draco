//----------------------------------*-C++-*----------------------------------//
// tstOpacity_pt.cc
// Thomas M. Evans
// Tue Apr 27 12:08:25 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Instantiations for tstOpacity
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "mc/RZWedge_Mesh.hh"
#include "../Opacity_Builder.t.hh"
#include "../Opacity.t.hh"
#include "../Mat_State.t.hh"
#include <iostream>

namespace rtt_imc
{

using rtt_mc::RZWedge_Mesh;
using rtt_mc::OS_Mesh;
using std::ostream;

// make instantiations of the objects that we want
template class Opacity_Builder<OS_Mesh>;

template class Opacity<OS_Mesh>;
template ostream& operator<<(ostream &, const Opacity<OS_Mesh> &);

template class Mat_State<OS_Mesh>;
template ostream& operator<<(ostream &, const Mat_State<OS_Mesh> &);

template class Mat_State<RZWedge_Mesh>;

} // end of namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of tstOpacity_pt.cc
//---------------------------------------------------------------------------//
