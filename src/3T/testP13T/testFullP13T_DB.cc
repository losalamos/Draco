#include "3T/testP13T/testFullP13T_DB.hh"

//----------------------------------*-C++-*----------------------------------//
// testFullP13T_DB.cc
// Randy M. Roberts
// Fri Jun 12 11:29:03 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/testP13T/testFullP13T_DB.hh"

#include "nml/Group.hh"
#include "nml/Items.hh"

using namespace XTM;

void testFullP13T_DB::setup_namelist( NML_Group& g )
{
#include ".nml_P13T.cc"
}



//---------------------------------------------------------------------------//
//                              end of testFullP13T_DB.cc
//---------------------------------------------------------------------------//
