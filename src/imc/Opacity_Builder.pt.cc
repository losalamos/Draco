//----------------------------------*-C++-*----------------------------------//
// Opacity_Builder.pt
// Thomas M. Evans
// Tue Mar 10 14:48:45 1998
//---------------------------------------------------------------------------//
// @> Instantialyzer for template class Opacity_Builder
//---------------------------------------------------------------------------//

#include "imctest/OS_Mesh.hh"
#include "imctest/OS_Interface.hh"
#include "ds++/SP.hh"
#include "Opacity_Builder.cc"

IMCSPACE

template class Opacity_Builder<OS_Mesh>;

template Opacity_Builder<OS_Mesh>::Opacity_Builder(SP<OS_Interface>, 
						   SP<OS_Mesh>);

CSPACE

//---------------------------------------------------------------------------//
//                              end of Opacity_Builder.pt
//---------------------------------------------------------------------------//
