//----------------------------------*-C++-*----------------------------------//
// Output_edit.hh
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Define the source.
//---------------------------------------------------------------------------//

#ifndef __sn_Output_edit_hh__
#define __sn_Output_edit_hh__

#include "sn/precision.hh"
#include "sn/Array.hh"
#include "sn/Input_edit.hh"

class Output_edit
{

    public:

        // use the default constructor
        // Output_edit();

        // use the default destructor
        // ~Output_edit();

        void print_flux_sum( Input_edit &data, Array4D &flux );

        void print_leakages( REAL *lkgs_l );

        void print_timing( REAL cpu0, REAL cpu1, REAL wall0,
                           REAL wall1, Input_edit &data, int its );

        void print_flux_moments( Input_edit &data, Array4D &flux );

    private:

        REAL flux_sum;  // a spatial sum of the scalar flux

};

#endif                          // __sn_Output_edit_hh__

//---------------------------------------------------------------------------//
//                              end of Output_edit.hh
//---------------------------------------------------------------------------//

