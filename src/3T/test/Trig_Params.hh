//----------------------------------*-C++-*----------------------------------//
// Trig_Params.hh
// Scott A. Turner (based on Geoffrey Furnish's Quad_Params.hh)
// 11 December 1997
//---------------------------------------------------------------------------//
// @> Hold onto parameters for the trigonometric test case.
//---------------------------------------------------------------------------//

#ifndef __3T_test_Trig_Params_hh__
#define __3T_test_Trig_Params_hh__

class NML_Group;

//===========================================================================//
// class Trig_Params - Parameters for the trigonometric test case

// This class holds the parameters needed to parameterize the trigonometric
// test case.  They are read in via the namelist utility.
//===========================================================================//

class Trig_Params {
  protected:

#include ".nml_Trig_Params.hh"

  public:
    void setup_namelist( NML_Group& g );
};

#endif                          // __3T/test_Trig_Params_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/test/Trig_Params.hh
//---------------------------------------------------------------------------//
