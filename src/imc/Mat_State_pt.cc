//----------------------------------*-C++-*----------------------------------//
// Mat_State_pt.cc
// Thomas M. Evans
// Tue Mar 10 14:33:24 1998
//---------------------------------------------------------------------------//
// @> Instantialyzer for template class Mat_State
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "Mat_State.t.hh"
#include <iostream>

IMCSPACE

using std::ostream;

template class Mat_State<OS_Mesh>;

template ostream& operator<<(ostream &, const Mat_State<OS_Mesh> &);

CSPACE

//---------------------------------------------------------------------------//
//                              end of Mat_State_pt.cc
//---------------------------------------------------------------------------//
