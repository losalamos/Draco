//----------------------------------*-C++-*----------------------------------//
// Host_Manager_pt.cc
// Thomas M. Evans
// Mon Aug  3 10:57:02 1998
//---------------------------------------------------------------------------//
// @> Instantialyzer for template class Host_Manager
//---------------------------------------------------------------------------//

#include "OS_Mesh.hh"
#include "AMR_Builder.hh"
#include "AMR_Interface.hh"
#include "Host_Manager.t.hh"

IMCSPACE

template class Host_Manager<OS_Mesh, AMR_Builder, AMR_Interface>;

CSPACE

//---------------------------------------------------------------------------//
//                              end of Host_Manager_pt.cc
//---------------------------------------------------------------------------//
