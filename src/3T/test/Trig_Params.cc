//----------------------------------*-C++-*----------------------------------//
// Trig_Params.cc
// Scott A. Turner (based on Geoffrey Furnish's Quad_Params.cc)
// 11 December 1997
//---------------------------------------------------------------------------//
// @> Hold onto parameters for the trigonometric test case.
//---------------------------------------------------------------------------//

#include "3T/test/Trig_Params.hh"

#include "nml/Group.hh"
#include "nml/Items.hh"

void Trig_Params::setup_namelist( NML_Group& g )
{
#include ".nml_Trig_Params.cc"
}

//---------------------------------------------------------------------------//
//                              end of Trig_Params.cc
//---------------------------------------------------------------------------//
