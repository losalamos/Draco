//----------------------------------*-C++-*----------------------------------//
// sweep3d.hh
// Scott Turner
// 17 March 1998
//---------------------------------------------------------------------------//
// @> Performs structured sweep.
//---------------------------------------------------------------------------//

//
// This routine performs an ordered sweep. The order is as follows:
// along j, as k is incremented, followed by incrementing i, then mu, then
// angle, then eta, then tsi.
//

#ifndef __sn_sweep3d_hh__
#define __sn_sweep3d_hh__

#include "sn/precision.hh"
#include "sn/array.hh"

class sweep3d
{

  public:

    // use the default constructor
    // sweep3d();

    // use the default destructor
    // ~sweep3d();

    void do_sweep(       int it,    int jt,  int kt,     int mm,  int nm,
                         int ibl,            int ibb,             int ibfr,
                         int ifxg,          REAL dx,             REAL dy,
                        REAL dz,            REAL *lkgs_l,  const REAL *w,
                  const REAL *mu,     const REAL *eta,     const REAL *tsi,
                  const REAL *wmu,    const REAL *weta,    const REAL *wtsi,
               const Array3D &p,   const Array3D &ct,         Array4D &src,
                     Array4D &flux                                         );

  private:

};

#endif                          // __sn_sweep3d_hh__

//---------------------------------------------------------------------------//
//                              end of sweep3d.hh
//---------------------------------------------------------------------------//

