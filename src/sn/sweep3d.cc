//----------------------------------*-C++-*----------------------------------//
// sweep3d.cc
// Scott Turner
// 19 February 1998
//---------------------------------------------------------------------------//
// @> Performs structured sweep.
//---------------------------------------------------------------------------//

//
// This routine performs an ordered sweep. The order is as follows:
// along j, as k is incremented, followed by incrementing i, then mu, then
// angle, then eta, then tsi.
//

#include <math.h>
#include <iostream.h>
#include <stdlib.h>

#include "sn/precision.hh"
#include "sn/protos.hh"

void sweep3d (       int it,    int jt,  int kt,     int mm,  int nm,
                     int ibl,            int ibb,             int ibfr,
                     int ifxg,          REAL dx,             REAL dy,
                    REAL dz,            REAL *lkgs_l,  const REAL *w,
              const REAL *mu,     const REAL *eta,     const REAL *tsi,
              const REAL *wmu,    const REAL *weta,    const REAL *wtsi,
           const Array3D &p,   const Array3D &ct,         Array4D &src,
                 Array4D &flux                                         )
{

#include "sn/array.hh"

  int itmm;     // number of x-direction cells * number of angles
  int i_cm;     // used in octant sweep order setup
  int mz;       // used in octant sweep order setup
  int ilim;     // used in octant sweep order setup
  int ishift;   // used in octant sweep order setup
  int jshift;   // used in octant sweep order setup
  int kshift;   // used in octant sweep order setup
  int jlow;     // used in octant sweep order setup
  int jhigh;    // used in octant sweep order setup
  int klow;     // used in octant sweep order setup
  int khigh;    // used in octant sweep order setup
  int ioff;     // used in octant sweep order setup
  int nk1;      // used to set first index of phi when building
  int nk2;      //   the angular source and the flux moments
  int iop;      // loop index counter for sweep over pair of quadrants
  int k1;       // index for looping over tsi
  int ih;       // index for looping over eta
  int iqp;      // quadrant number(iq) + 1
  int iz;       // first index of phi during balance equation
  int i;        // loop variable for cells in the x-direction
  int j;        // loop variable for cells in the y-direction
  int k;        // loop variable for cells in the z-direction
  int m;        // loop variable for quadrature points (angles) per quadrant
  int n;        // loop variable for flux and source angular moments
  int iq;       // loop variable for the number of quadrants, may be calculated,
                //   rather than incremented in some cases
  int jbdim;    // dimension parameter for boundary arrays
  int kbdim;    // dimension parameter for boundary arrays
  int maxop;    // loop index for quadrant sweep

  REAL c1;                    // used to hold a value calculated outside of a
  REAL c2;                    //   loop, for use inside the loop (efficiency)
  REAL ql;                    // a temp variable to hold calculated results
  REAL dl;                    // a temp variable to hold calculated results
  REAL sih;                   // temp holder for phii during set-to-zero fixup
  REAL siv;                   // temp holder for phij during set-to-zero fixup
  REAL sif;                   // temp holder for phik during set-to-zero fixup

  Array2D phii(jt,kt);        // temp angular flux for the left face of a cell
  Array3D phi(jt,kt,it*mm*2); // cell centered angular flux
  Array3D fh_i(jt,kt,2);      // left and right leakages (horizontal,x-axis)
  Array3D fv_i(kt,it,2);      // bottom and top leakages (vertical  ,y-axis)
  Array3D fz_i(jt,it,2);      // front and back leakages (in-out    ,z-axis)
  Array4D dlinv(jt,kt,it,mm); // a factor in the balance equation, which is
                              //   pre-calculated outside of the loop over
                              //   quadrants, for efficiency

// Allocate memory for local arrays.

  REAL *const phij = new REAL [kt]; // temp angular flux for top  face of a cell
  REAL *const phik = new REAL [jt]; // temp angular flux for back face of a cell
  REAL *const hi   = new REAL [it]; // 2.0 / dx (a commonly used factor)
  REAL *const hj   = new REAL [jt]; // 2.0 / dy (a commonly used factor)
  REAL *const hk   = new REAL [kt]; // 2.0 / dz (a commonly used factor)

// Set boundary array dimensions based on boundary condition.
 
  if (ibb == 1)
    jbdim = 2*mm*it;                          // reflective boundary
  else
    jbdim = 1;                                // vacuum boundary
 
  if (ibfr == 1)
    kbdim = 2*mm*it;                          // reflective boundary
  else
    kbdim = 1;                                // vacuum boundary

// Initialize arrays. bsavv and bsavz hold the edge fluxes for reflective bdry
// conditions. bsavv is 2d because the 2nd pair of octants uses phij data from
// the 1st, and the 4th pair uses data from the 3rd, which is saved over the
// 1st. bsavz is 3d because it depends on phik data from the 1st pair while
// doing a 3rd pair calc, and from the 2nd pair during a 4th pair calc, thus
// overwriting is not possible.

  Array2D bsavv(kt,jbdim);
  Array3D bsavz(jt,kbdim,2);

  if ( ibb == 1 )  bsavv.Array2D_reinit(0.0);
  if ( ibfr == 1 ) bsavz.Array3D_reinit(0.0);

  fh_i.Array3D_reinit(0.0);
  fv_i.Array3D_reinit(0.0);
  fz_i.Array3D_reinit(0.0);

// Pre-calculate commonly used factors outside of the loops for efficiency.

  itmm = it * mm;

  for ( i=0 ; i < it ; i++ )
    hi[i] = 2.0 / dx;
  for ( j=0 ; j < jt ; j++ )
    hj[j] = 2.0 / dy;
  for ( k=0 ; k < kt ; k++ )
    hk[k] = 2.0 / dz;

  for ( m=0 ; m < mm ; m++ ) {
  for ( i=0 ; i < it ; i++ ) {
  for ( k=0 ; k < kt ; k++ ) {
  for ( j=0 ; j < jt ; j++ )
    dlinv(j,k,i,m) = 1.0 / ( ct(j,k,i)  +  mu[m] * hi[i]
                                        + eta[m] * hj[j]
                                        + tsi[m] * hk[k] );
  } } }

//******************************************************************************
// Begin the loop over quadrant pairs, where k1 is a loop over tsi, ih is a loop
// over eta, and mu is reversed within each ordered (thus creating a pair of
// quadrants).
//******************************************************************************

  for ( k1=0 ; k1 < 2 ; k1++ ) {
  for ( ih=0 ; ih < 2 ; ih++ ) {

    iq  = 2*ih + 4*k1;
    iqp = iq + 1;

// Set up the order in which octants will be swept. The order allows parallel
// sweeps to be ordered across (j,k) for stacks of i-values. This
// allows us to match the old vectorized sweepers' logic and order which was
// optimized for the fact that reflective boundaries typically occured at the
// left(i), bottom(j), and front(k) faces of most problems. Beginning sweeps
// toward reflective boundaries, followed by sweeps away from them reduces the
// amount of angular flux information that must be stored and give you explicit
// boundaries rather than implicit ones. jlow, jhigh, klow, and khigh
// specify the spatial sweep order within an i-plane (a plane normal to x-axis).
// During the loop over iop=1,maxop where maxop=2*mm*it, mu alternates between
// -1 and 1 rather than doing all mm angles for mu=-1 then all mm angles for
// mu=1. This allows for reuse of the i-direction's angular flux array (phii)
// at the reflective boundary which may be completely contained in cache when
// mu changes sign and the calculation proceeds back along the same angle m.
// More important than the cache affect though, is the avoidance of having to
// store phii for each angle at the reflective boundary.
//
// The order for the octant sweeps is:
//
//  iq   mu   eta   tsi
//   1   -1    -1    -1
//   2   +1    -1    -1
//   3   -1    +1    -1
//   4   +1    +1    -1
//   5   -1    -1    +1
//   6   +1    -1    +1
//   7   -1    +1    +1
//   8   +1    +1    +1 

    i_cm   = it + 1;
    ishift = -1;
    ioff   = 0;
    mz     = 0;
    ilim   = 1;

    if ( ih == 0 )
    {
      jshift = -1;
      jlow   = jt-1;
      jhigh  = 0;
    }
    else
    {
      jshift = +1;
      jlow   = 0;
      jhigh  = jt-1;
    }

    if ( k1 == 0 )
    {
      kshift = -1;
      klow   = kt-1;
      khigh  = 0;
    }
    else
    {
      kshift = +1;
      klow   = 0;
      khigh  = kt-1;
    }

// Build the angular source in phi. phi must first be re-initialized to zero.

    phi.Array3D_reinit(0.0);

    for ( n=0 ; n < nm ; n++ ) {
    for ( m=0 ; m < mm ; m++ ) {

      nk1 = m * it;
      nk2 = nk1 + mm * it;
      c1 = p(n,m,iq);
      c2 = p(n,m,iqp);

      for ( i=0 ; i < it ; i++ ) {
      for ( k=0 ; k < kt ; k++ ) {
      for ( j=0 ; j < jt ; j++ ) {
        phi(j,k,i+nk1) += src(j,k,i,n) * c1;
        phi(j,k,i+nk2) += src(j,k,i,n) * c2;
      } } }

    } }

//******************************************************************************
// Begin the sweep over a pair of quadrants.
//******************************************************************************

    maxop  = 2*mm*it;
    for ( iop=0 ; iop < maxop ; iop++ ) {

// Update sweep indices and control variables.

      i_cm += ishift;
      i = i_cm - 1;

// Due to the logic flow of this loop, we need to test for i > it-1, which in
// C is equivalent to the it+1 cell (outside of physical system). If this occurs
// then skip ahead to the sweep reversal.

      if ( i <= it-1 )
      {
        iz = i + mz*it + itmm*ioff;

// Initialize edge fluxes, impose boundary conditions. Leakage is calculated
// for reflective boundaries here to cancel out (via opposite sign) a later
// calculation that is done at all boundaries. It is assumed here, that the
// right, top, and back boundaries are vacuum boundaries. If any of the left,
// bottom, or front boundaries are also vacuum, then their corresponding edge
// flux arrays are set to zero at the time of sweep reversal. If the bottom or
// front boundaries are reflective than phij or phik must be loaded from bsavv
// or basavz at the time of sweep reversal.

        if ( ishift == -1 && i == it-1 )
 
          phii.Array2D_reinit(0.0);
 
        else if ( ishift == +1 && i == 0 && ibl == 0 )

          phii.Array2D_reinit(0.0);

        if ( ih == 0 )                            // corners 1,2,5,6
        {
          for ( k=0 ; k < kt ; k++ )
            phij[k] = 0.0;
        }
        else                                      // corners 3,4,7,8
        {
          if ( ibb == 0 )
          {
            for ( k=0 ; k < kt ; k++ )
              phij[k] = 0.0;
          }
          else if ( ibb == 1 )
          {
            for ( k=0 ; k < kt ; k++ ) {
              phij[k] = bsavv(k,iop);
            }
          }
        }

        if ( k1 == 0 )                            // corners 1,2,3,4
        {
          for ( j=0 ; j < jt ; j++ )
            phik[j] = 0.0;
        }
        else                                      // corners 5,6,7,8
        {
          if ( ibfr == 0 )
          {
            for ( j=0 ; j < jt ; j++ )
              phik[j] = 0.0;
          }
          else if ( ibfr == 1 )
          {
            for ( j=0 ; j < jt ; j++ ) {
              phik[j] = bsavz(j,iop,ih);
            }
          }
        }

        if ( ifxg == 0 )
        {

// Balance equation with diamond difference and no flux fixup. phii, phij, phik
// are only needed from one i-plane to the next and thus can be continuously
// overwritten. phij and phik are saved/restored to/from bsavv and bsavz at
// boundaries. phii is immediately reused via sweep reversal and thus does not
// require a storage array.

          for ( k=klow ; k != (khigh+kshift) ; k += kshift ) {
          for ( j=jlow ; j != (jhigh+jshift) ; j += jshift ) {

            ql = phi(j,k,iz) +
                 mu [mz] * phii(j,k) * hi[i] +
                 eta[mz] * phij[k]   * hj[j] +
                 tsi[mz] * phik[j]   * hk[k];

            phi(j,k,iz) = ql * dlinv(j,k,i,mz);

            phii(j,k) = 2.0 * phi(j,k,iz) - phii(j,k);
            phij[k]   = 2.0 * phi(j,k,iz) - phij[k];
            phik[j]   = 2.0 * phi(j,k,iz) - phik[j];

          } }
        }
        else
        {

// Balance equation w/ diamond difference and set to zero flux fixup.

          for ( k=klow ; k != (khigh+kshift) ; k += kshift ) {
          for ( j=jlow ; j != (jhigh+jshift) ; j += jshift ) {

            ql = phi(j,k,iz) +
                 mu [mz] * phii(j,k) * hi[i] +
                 eta[mz] * phij[k]   * hj[j] +
                 tsi[mz] * phik[j]   * hk[k];
            dl = ct(j,k,i) + 
                 mu [mz] * hi[i] +
                 eta[mz] * hj[j] +
                 tsi[mz] * hk[k];

            phi(j,k,iz) = ql * dlinv(j,k,i,mz);

            sih = 2.0 * phi(j,k,iz) - phii(j,k);
            siv = 2.0 * phi(j,k,iz) - phij[k];
            sif = 2.0 * phi(j,k,iz) - phik[j];

// Test for negative fluxes, if found, adjust temporary value holders,
// ql and dl, to account for setting the appropriate negative flux to
// zero. After setting a flux to zero and adjusting other fluxes to
// account for it, then go back to make sure the change did not cause
// other fluxes to become negative.

  neg_flux:

            if ( sih < 0.0 )
            {
              ql -= 0.5 * mu[mz] * phii(j,k) * hi[i];
              dl -= mu[mz] * hi[i];

              phi(j,k,iz) = ql / dl;

              sih = 0.0;
              if ( siv != 0.0 ) siv = 2.0 * phi(j,k,iz) - phij[k];
              if ( sif != 0.0 ) sif = 2.0 * phi(j,k,iz) - phik[j];
            }

            if ( siv < 0.0 )
            {
              ql -= 0.5 * eta[mz] * phij[k] * hj[j];
              dl -= eta[mz] * hj[j];

              phi(j,k,iz) = ql / dl;

              siv = 0.0;
              if ( sif != 0.0 ) sif = 2.0 * phi(j,k,iz) - phik[j];
              if ( sih != 0.0 ) sih = 2.0 * phi(j,k,iz) - phii(j,k);

              goto neg_flux;
            }

            if ( sif < 0.0 )
            {
              ql -= 0.5 * tsi[mz] * phik[j] * hk[k];
              dl -= tsi[mz] * hk[k];

              phi(j,k,iz) = ql / dl;

              sif = 0.0;
              if ( sih != 0.0 ) sih = 2.0 * phi(j,k,iz) - phii(j,k);
              if ( siv != 0.0 ) siv = 2.0 * phi(j,k,iz) - phij[k];

              goto neg_flux;
            }

            phii(j,k) = sih;
            phij[k]   = siv;
            phik[j]   = sif;

          } }  // end of flux fixup j,k loops

        }  // end of flux fixup conditional

// Save the j and k edge fluxes at reflective boundaries. The i-edge fluxes do
// not need to be save because they are about to be used in the sweep reversal.

        if ( ibb  == 1 && ih == 0 )
        {
          for ( k=klow ; k != (khigh+kshift) ; k += kshift )
            bsavv(k,iop) = phij[k];
        }

        if ( ibfr == 1 && k1 == 0 )
        {
          for ( j=jlow ; j != (jhigh+jshift) ; j += jshift )
            bsavz(j,iop,ih) = phik[j];
        }

// Compute leakages.

        if ( i+ioff == 0 && ibl == 0 )
        {
          c1 = ishift * wmu[mz];
          for ( k=0 ; k < kt ; k++ ) {
          for ( j=0 ; j < jt ; j++ )
            fh_i(j,k,0) += c1 * phii(j,k);
          }
        }

        if ( i+ioff == it ) 
        {
          c1 = ishift * wmu[mz];
          for ( k=0 ; k < kt ; k++ ) {
          for ( j=0 ; j < jt ; j++ )
            fh_i(j,k,1) += c1 * phii(j,k);
          }
        }

        if ( ih == 1 || ibb == 0 )
        {
          c1 = jshift * weta[mz];
          for ( k=0 ; k < kt ; k++ )
            fv_i(k,i,ih) += c1 * phij[k];
        }

        if ( k1 == 1 || ibfr == 0 )
        {
          c1 = kshift * wtsi[mz];
          for ( j=0 ; j < jt ; j++ )
            fz_i(j,i,k1) += c1 * phik[j];
        }

      }  // end of (i<it-1) conditional

//******************************************************************************
// Sweep reversal for mu direction cosine.
//******************************************************************************

      if ( i_cm == ilim )
      {
        if ( ishift == -1 )
        {
          i_cm   = 0;
          ishift = +1;
          ioff   = +1;
          ilim   = it;
        }
        else if ( mz < mm-1 )
        {
          i_cm   = it + 1;
          ishift = -1;
          ioff   = 0;
          ilim   = 1;
          mz    += 1;
        }
      }

    }  // end of iop loop

//******************************************************************************
// End of the sweep over a pair of quadrants.
//******************************************************************************

// Compute the flux moments.

    for ( m=0 ; m < mm ; m++ ) {

      nk1 = m * it;
      nk2 = nk1 + itmm;

      for ( i=0 ; i < it ; i++ ) {
      for ( k=0 ; k < kt ; k++ ) {
      for ( j=0 ; j < jt ; j++ )
        flux(j,k,i,0) += w[m]*(phi(j,k,i+nk1) + phi(j,k,i+nk2));
      } }
    }

    for ( n=1 ; n < nm ; n++ ) {
    for ( m=0 ; m < mm ; m++ ) {

      nk1 = m * it;
      nk2 = nk1 + itmm;
      c1  = w[m] * p(n,m,iq);
      c2  = w[m] * p(n,m,iqp);

      for ( i=0 ; i < it ; i++ ) {
      for ( k=0 ; k < kt ; k++ ) {
      for ( j=0 ; j < jt ; j++ )
        flux(j,k,i,n) += c1*phi(j,k,i+nk1) + c2*phi(j,k,i+nk2);
      } }

    } }

  } }  // end of k1 and ih loops

//******************************************************************************
// End of the loop over quadrant pairs.
//******************************************************************************

// Calculate leakages for the balance table.

  lkgs_l[0] = 0.0;  // left+right leakage (also known as horizontal or x-axis)
  lkgs_l[1] = 0.0;  // right      leakage
  lkgs_l[2] = 0.0;  // bottom+top leakage (also known as vertical or y-axis)
  lkgs_l[3] = 0.0;  // top        leakage
  lkgs_l[4] = 0.0;  // front+back leakage (also known as in-out or z-axis)
  lkgs_l[5] = 0.0;  // back       leakage

  for ( k=0 ; k < kt ; k++ ) {
  for ( j=0 ; j < jt ; j++ ) {

    fh_i(j,k,0) *= 4.0 / (hj[j] * hk[k]);
    fh_i(j,k,1) *= 4.0 / (hj[j] * hk[k]);
    lkgs_l[0] -= fh_i(j,k,0);
    lkgs_l[1] += fh_i(j,k,1);

  } }

  for ( i=0 ; i < it ; i++ ) {
  for ( k=0 ; k < kt ; k++ ) {

    fv_i(k,i,0) *= 4.0 / (hi[i] * hk[k]);
    fv_i(k,i,1) *= 4.0 / (hi[i] * hk[k]);
    lkgs_l[2] -= fv_i(k,i,0);
    lkgs_l[3] += fv_i(k,i,1);

  } }

  for ( i=0 ; i < it ; i++ ) {
  for ( j=0 ; j < jt ; j++ ) {

    fz_i(j,i,0) *= 4.0 / (hi[i] * hj[j]);
    fz_i(j,i,1) *= 4.0 / (hi[i] * hj[j]);
    lkgs_l[4] -= fz_i(j,i,0);
    lkgs_l[5] += fz_i(j,i,1);

  } }

  lkgs_l[0] += lkgs_l[1];
  lkgs_l[2] += lkgs_l[3];
  lkgs_l[4] += lkgs_l[5];

  delete [] phij;
  delete [] phik;
  delete [] hi  ;
  delete [] hj  ;
  delete [] hk  ;

  return;
}

//---------------------------------------------------------------------------//
//                              end of sweep3d.cc
//---------------------------------------------------------------------------//

