//----------------------------------*-C++-*----------------------------------//
// IMC_Manager_pt.cc
// Thomas M. Evans
// Mon Jun 15 17:03:01 1998
//---------------------------------------------------------------------------//
// @> Instantialyzer for template class IMC_Manager
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "OS_Builder.hh"
#include "OS_Interface.hh"
#include "IMC_Manager.t.hh"

IMCSPACE

template class IMC_Manager<OS_Mesh, OS_Builder, OS_Interface>;

CSPACE

//---------------------------------------------------------------------------//
//                              end of IMCTEST_Man_pt.cc
//---------------------------------------------------------------------------//
