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
#include "sn/Input_edit.hh"

#include "ds++/Mat.hh"
using dsxx::Mat3;
using dsxx::Mat4;

class Error
{

    public:

        // the constructor is used to initialize private data

        Error( Input_edit &data );

        // use the default destructor
        // ~Error();

        void error_init( Input_edit &data, Mat4<REAL> &flux );
        void error_calc( Input_edit &data, Mat4<REAL> &flux );
        void error_print( int its );
        bool error_test( Input_edit &data, int its );

    private:

        Mat3<REAL> pflux;  // both the previous flux and the relative change in
                           //   the flux after each iteration
        REAL err;          // max change in flux from one iteration to the next

};

#endif                          // __sn_Error_hh__

//---------------------------------------------------------------------------//
//                              end of Error.hh
//---------------------------------------------------------------------------//

