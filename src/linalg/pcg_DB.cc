//----------------------------------*-C++-*----------------------------------//
// pcg_DB.cc
// Dave Nystrom
// Fri May  9 09:45:08 1997
//---------------------------------------------------------------------------//
// @> pcg descriptor block.
//---------------------------------------------------------------------------//

#include "pcg_DB.hh"

//---------------------------------------------------------------------------//
// Attach the parameters to a namelist group.
//---------------------------------------------------------------------------//

void pcg_DB::setup_namelist( NML_Group& g )
{
    dsxx::String fldeqn = name;

#include ".nml_pcg.cc"
}

//---------------------------------------------------------------------------//
//                              end of pcg_DB.cc
//---------------------------------------------------------------------------//
