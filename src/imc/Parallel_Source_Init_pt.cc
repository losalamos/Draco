//----------------------------------*-C++-*----------------------------------//
// Parallel_Source_Init_pt.cc
// Thomas M. Evans
// Fri Aug 21 17:22:58 1998
//---------------------------------------------------------------------------//
// @> Instantialyzer for template class Parallel_Source_Init
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "AMR_Interface.hh"
#include "ds++/SP.hh"
#include "Parallel_Source_Init.t.hh"

IMCSPACE

using dsxx::SP;

template class Parallel_Source_Init<OS_Mesh>;

template 
Parallel_Source_Init<OS_Mesh>::Parallel_Source_Init(SP<AMR_Interface>, 
						    SP<OS_Mesh>);

CSPACE

//---------------------------------------------------------------------------//
//                              end of Parallel_Source_Init_pt.cc
//---------------------------------------------------------------------------//
