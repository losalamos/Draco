//----------------------------------*-C++-*----------------------------------//
// Error.hh
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Error convergence.
//---------------------------------------------------------------------------//

#ifndef __sn_Error_hh__
#define __sn_Error_hh__

#include "sn/precision.hh"
#include "sn/Array.hh"
#include "sn/Input_edit.hh"

class Error
{

    public:

        // the constructor is used to initialize private data

        Error( Input_edit &data );

        // use the default destructor
        // ~Error();

        void error_init( Input_edit &data, Array4D &flux );
        void error_calc( Input_edit &data, Array4D &flux );
        void error_print( int its );
        bool error_test( Input_edit &data, int its );

    private:

        Array3D pflux;  // both the previous flux and the relative change in the
                        //   flux after each iteration
        REAL err;       // max change in the flux from one iteration to the next

};

#endif                          // __sn_Error_hh__

//---------------------------------------------------------------------------//
//                              end of Error.hh
//---------------------------------------------------------------------------//

