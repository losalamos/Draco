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

IMCSPACE

using std::ostream;

template class Opacity<OS_Mesh>;

template ostream& operator<<(ostream&, const Opacity<OS_Mesh> &);

CSPACE

//---------------------------------------------------------------------------//
//                              end of Opacity_pt.cc
//---------------------------------------------------------------------------//
