//----------------------------------*-C++-*----------------------------------//
// Quad_Params.hh
// Geoffrey M. Furnish
// Thu Nov 20 11:11:27 1997
//---------------------------------------------------------------------------//
// @> Hold onto parameters for the quadratic test case.
//---------------------------------------------------------------------------//

#ifndef __3T_test_Quad_Params_hh__
#define __3T_test_Quad_Params_hh__

class NML_Group;

//===========================================================================//
// class Quad_Params - Parameters for the quadratic test case

// This class holds the parameters needed to parameterize the quadratic test
// case.  They are read in via the namelist utility.
//===========================================================================//

class Quad_Params {
  protected:

#include ".nml_Quad_Params.hh"

  public:
    void setup_namelist( NML_Group& g );
};

#endif                          // __3T/test_Quad_Params_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/test/Quad_Params.hh
//---------------------------------------------------------------------------//
