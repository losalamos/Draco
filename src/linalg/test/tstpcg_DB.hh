//----------------------------------*-C++-*----------------------------------//
// tstpcg_DB.hh
// Dave Nystrom
// Mon May 12 10:52:19 1997
//---------------------------------------------------------------------------//
// @> tstpcg descriptor block.
//---------------------------------------------------------------------------//

#ifndef __linalg_test_tstpcg_DB_hh__
#define __linalg_test_tstpcg_DB_hh__

#include "nml/Group.hh"
#include "nml/Items.hh"

//===========================================================================//
// class tstpcg_DB - tstpcg descriptor block

// Holds namelist parameters for setting up and controlling tstpcg.
//===========================================================================//

class tstpcg_DB {
  public:

// Pick up the namelist variables.

#include ".nml_tstpcg.hh"

// Methods.

  public:
    void setup_namelist( NML_Group& g );

};

#endif                          // __linalg_test_tstpcg_DB_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/test/tstpcg_DB.hh
//---------------------------------------------------------------------------//
