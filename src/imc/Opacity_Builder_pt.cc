//----------------------------------*-C++-*----------------------------------//
// Opacity_Builder_pt.cc
// Thomas M. Evans
// Tue Mar 10 14:48:45 1998
//---------------------------------------------------------------------------//
// @> Instantialyzer for template class Opacity_Builder
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "OS_Interface.hh"
#include "AMR_Interface.hh"
#include "ds++/SP.hh"
#include "Opacity_Builder.t.hh"

IMCSPACE

using dsxx::SP;

template class Opacity_Builder<OS_Mesh>;

template Opacity_Builder<OS_Mesh>::Opacity_Builder(SP<OS_Interface>);

template Opacity_Builder<OS_Mesh>::Opacity_Builder(SP<AMR_Interface>);

CSPACE

//---------------------------------------------------------------------------//
//                              end of Opacity_Builder_pt.cc
//---------------------------------------------------------------------------//
