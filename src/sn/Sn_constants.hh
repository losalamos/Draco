//----------------------------------*-C++-*----------------------------------//
// Sn_constants.hh
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Define the Sn constants.
//---------------------------------------------------------------------------//

#ifndef __sn_Sn_constants_hh__
#define __sn_Sn_constants_hh__

#include "sn/precision.hh"
#include "sn/Array.hh"
#include "sn/Input_edit.hh"

class Sn_constants
{

    public:

        // the ctor is used to set quadrature and spherical harmonics constants

        Sn_constants( Input_edit &data );

        // define a destructor

        ~Sn_constants();

        REAL w(    int m ) const { return w_p[   m]; }
        REAL mu(   int m ) const { return mu_p[  m]; }
        REAL eta(  int m ) const { return eta_p[ m]; }
        REAL tsi(  int m ) const { return tsi_p[ m]; }
        REAL wmu(  int m ) const { return wmu_p[ m]; }
        REAL weta( int m ) const { return weta_p[m]; }
        REAL wtsi( int m ) const { return wtsi_p[m]; }

        REAL p( int n, int m, int iq ) const { return p_p(n,m,iq); }

    private:

        REAL *w_p;
        REAL *mu_p;
        REAL *eta_p;
        REAL *tsi_p;
        REAL *wmu_p;
        REAL *weta_p;
        REAL *wtsi_p;

        Array3D p_p;  // fixed source

};

#endif                          // __sn_Sn_constants_hh__

//---------------------------------------------------------------------------//
//                              end of Sn_constants.hh
//---------------------------------------------------------------------------//

