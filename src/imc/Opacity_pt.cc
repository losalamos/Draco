//----------------------------------*-C++-*----------------------------------//
// Opacity_pt.cc
// Thomas M. Evans
// Tue Mar 10 14:47:45 1998
//---------------------------------------------------------------------------//
// @> Instantialyzer for template class Opacity
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "Opacity.t.hh"
#include <iostream>

namespace rtt_imc 
{

using std::ostream;

template class Opacity<OS_Mesh>;

template ostream& operator<<(ostream&, const Opacity<OS_Mesh> &);

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Opacity_pt.cc
//---------------------------------------------------------------------------//
