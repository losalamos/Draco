//----------------------------------*-C++-*----------------------------------//
// Source_Init_pt.cc
// Thomas M. Evans
// Wed May  6 15:05:32 1998
//---------------------------------------------------------------------------//
// @> Instantialyzer for template class Source_Init
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "OS_Interface.hh"
#include "AMR_Interface.hh"
#include "ds++/SP.hh"
#include "Source_Init.t.hh"

IMCSPACE

using dsxx::SP;

template class Source_Init<OS_Mesh>;

template Source_Init<OS_Mesh>::Source_Init(SP<OS_Interface>, SP<OS_Mesh>);

template Source_Init<OS_Mesh>::Source_Init(SP<AMR_Interface>, SP<OS_Mesh>);

CSPACE

//---------------------------------------------------------------------------//
//                              end of Source_Init_pt.cc
//---------------------------------------------------------------------------//
