//----------------------------------*-C++-*----------------------------------//
// read_input.hh
// Scott Turner
// 17 March 1998
//---------------------------------------------------------------------------//
// @> Read and test input.
//---------------------------------------------------------------------------//

#ifndef __sn_read_input_hh__
#define __sn_read_input_hh__

#include "sn/precision.hh"

class read_input
{

  public:

    // use the default destructor
    // read_input();

    // use the default destructor
    // ~read_input();

    void read_data();

    void get_all_data( int &it_l, int &jt_l, int &kt_l, int &mm_l,
                       int &isct_l, int &ibl_l, int &ibb_l,
                       int &ibfr_l, int &iprint_l, int &ifxg_l,
                       REAL &dx_l, REAL &dy_l, REAL &dz_l,
                       REAL &epsi_l );

    void get_basic_data( int &it_l, int &jt_l, int &kt_l, int &mm_l,
                         int &isct_l );

  private:

    int it;       // total number of mesh cells in the x-direction
    int jt;       // total number of mesh cells in the y-direction
    int kt;       // total number of mesh cells in the z-direction
    int mm;       // number of   quadrature points (angles) per quadrant
    int isct;     // legendre order of scattering
    int ibl;      // left boundary condition 0/1 = vacuum/reflective
    int ibb;      // bottom boundary condition 0/1 = vacuum/reflective
    int ibfr;     // front boundary condition 0/1 = vacuum/reflective
    int iprint;   // 0/1 = no/yes print fluxes
    int ifxg;     // 0/1 = no/yes set to zero flux fixup

    REAL dx;      // mesh size in the x-direction
    REAL dy;      // mesh size in the y-direction
    REAL dz;      // mesh size in the z-direction
    REAL epsi;    // convergence precision or,
                  //   if negative, then the number of iterations to do

};

#endif                          // __sn_read_input_hh__

//---------------------------------------------------------------------------//
//                              end of read_input.hh
//---------------------------------------------------------------------------//

