//----------------------------------*-C++-*----------------------------------//
// protos.hh
// Scott Turner
// 19 February 1998
//---------------------------------------------------------------------------//
// @> Function Prototypes
//---------------------------------------------------------------------------//

#ifndef __sn_protos_hh__
#define __sn_protos_hh__

//===========================================================================//
// Here are the function prototypes for sn++.
//===========================================================================//

#include "sn/precision.hh"
#include "sn/array.hh"

// Solver kernel

extern void sweep3d (       int it,    int jt,  int kt,     int mm,  int nm,
                            int ibl,            int ibb,             int ibfr,
                            int ifxg,          REAL dx,             REAL dy,
                           REAL dz,            REAL *lkgs_l,  const REAL *w,
                     const REAL *mu,     const REAL *eta,     const REAL *tsi,
                     const REAL *wmu,    const REAL *weta,    const REAL *wtsi,
                  const Array3D &p,   const Array3D &ct,         Array4D &src,
                        Array4D &flux                                         );

// CPU and Wallclock timing routine

extern void timer( REAL &cpu, REAL &wall );

#endif                          // __sn_protos_hh__

//---------------------------------------------------------------------------//
//                              end of protos.hh
//---------------------------------------------------------------------------//

