//----------------------------------*-C++-*----------------------------------//
// Quad_Params.cc
// Geoffrey M. Furnish
// Thu Nov 20 11:11:28 1997
//---------------------------------------------------------------------------//
// @> Hold onto parameters for the quadratic test case.
//---------------------------------------------------------------------------//

#include "3T/test/Quad_Params.hh"

#include "nml/Group.hh"
#include "nml/Items.hh"

void Quad_Params::setup_namelist( NML_Group& g )
{
#include ".nml_Quad_Params.cc"
}

//---------------------------------------------------------------------------//
//                              end of Quad_Params.cc
//---------------------------------------------------------------------------//
