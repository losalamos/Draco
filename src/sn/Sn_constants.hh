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
#include "sn/Input_edit.hh"

#include "ds++/Mat.hh"
using dsxx::Mat1;
using dsxx::Mat3;

class Sn_constants
{

    public:

        // the ctor is used to set quadrature and spherical harmonics constants

        Sn_constants( Input_edit &data );

        // use the default destructor
        // ~Sn_constants();

        REAL w(    int m ) const { return w_p(   m); }
        REAL mu(   int m ) const { return mu_p(  m); }
        REAL eta(  int m ) const { return eta_p( m); }
        REAL tsi(  int m ) const { return tsi_p( m); }
        REAL wmu(  int m ) const { return wmu_p( m); }
        REAL weta( int m ) const { return weta_p(m); }
        REAL wtsi( int m ) const { return wtsi_p(m); }

        REAL p( int n, int m, int iq ) const { return p_p(n,m,iq); }

    private:

        Mat1<REAL> w_p;
        Mat1<REAL> mu_p;
        Mat1<REAL> eta_p;
        Mat1<REAL> tsi_p;
        Mat1<REAL> wmu_p;
        Mat1<REAL> weta_p;
        Mat1<REAL> wtsi_p;

        Mat3<REAL> p_p;  // fixed source

};

#endif                          // __sn_Sn_constants_hh__

//---------------------------------------------------------------------------//
//                              end of Sn_constants.hh
//---------------------------------------------------------------------------//

