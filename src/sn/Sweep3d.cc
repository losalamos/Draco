//----------------------------------*-C++-*----------------------------------//
// Sweep3d.cc
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Performs structured sweep.
//---------------------------------------------------------------------------//

#include "sn/Sweep3d.hh"

void Sweep3d::do_sweep( Input_edit &data, Cross_section &xsec, Sn_constants &sn,
                        Pre_calcs &pre,   REAL *lkgs_l,        Array4D &src_mom,
                        Array4D &flux                                          )
{
    // This routine performs an ordered sweep. The order is as follows:
    // along j, as k is incremented, followed by incrementing i, then mu, then
    // angle, then eta, then tsi.

    // Initialize edge flux and boundary arrays.

    Array3D fh_i(data.jt(),data.kt(),2);  // left/right leakage (horiz,  x-axis)
    Array3D fv_i(data.kt(),data.it(),2);  // bottom/top leakage (vert ,  y-axis)
    Array3D fz_i(data.jt(),data.it(),2);  // front/back leakage (in-out, z-axis)
    Array2D bsavv(data.kt(),data.jbdim());    // storage of vertical boundaries
    Array3D bsavz(data.jt(),data.kbdim(),2);  // storage of in-out   boundaries

    //**************************************************************************
    // Begin the loop over quadrant pairs, where k1 is a loop over tsi, ih is a
    // loop over eta, and mu is reversed within each ordered (thus creating a
    // pair of quadrants).
    //**************************************************************************

    REAL *phij = new REAL [data.kt()]; // temp angular flux, top  cell face
    REAL *phik = new REAL [data.jt()]; // temp angular flux, back cell face

    Array2D phii(data.jt(),data.kt()); // temp angular flux, left cell face
    Array3D phi(data.jt(),data.kt(),data.it()*data.mm()*2); // cell centered
                                                            // angular flux
    for ( k1=0 ; k1 < 2 ; k1++ )
    for ( ih=0 ; ih < 2 ; ih++ )
    {
        iq  = 2 * ih + 4 * k1;
        iqp = iq + 1;

        // Set up the order in which octants will be swept.

        octant_ordering( data );

        // Build the angular source in phi.

        build_angular_source( data, sn, src_mom, phi );

        //**********************************************************************
        // Begin the sweep over a pair of quadrants.
        //**********************************************************************

        for ( iop=0 ; iop < data.maxop() ; iop++ )
        {
            // Update sweep indices and control variables.

            i_cm += ishift;
            isw = i_cm - 1;

            // Due to the logic flow of this loop, we need to test for
            // isw > it-1, which in C is equivalent to the it+1 cell (outside
            // of the physical system). If this occurs then skip ahead to the
            // sweep reversal.

            if ( isw <= data.it()-1 )
            {
                iz = isw + mz * data.it() + data.itmm() * ioff;

                // Initialize edge fluxes, impose boundary conditions.

                edge_and_boundary_set( data, phii, phij, phik, bsavv, bsavz );

                // Choose the balance equation to be solved.

                if ( data.ifxg() == 0 )
                {
                    balance_eqn_no_fixup( data, sn,   pre, phi,
                                          phii, phij, phik      );
                }
                else
                {
                    balance_eqn_with_fixup( data, sn,   pre,  xsec,
                                            phi,  phii, phij, phik  );
                }

                // Save the j and k edge fluxes at reflective boundaries.

                save_boundary( data, bsavv, bsavz, phij, phik );

                // Compute leakages.

                cell_boundary_leakage( data, sn,   fh_i, fv_i,
                                       fz_i, phii, phij, phik  );
            }

            //******************************************************************
            // Sweep reversal for mu direction cosine.
            //******************************************************************

            sweep_reversal( data );
        }

        //**********************************************************************
        // End of the sweep over a pair of quadrants.
        //**********************************************************************

        // Compute the flux moments.

        flux_moments( data, sn, flux, phi );
    }

    //**************************************************************************
    // End of the loop over quadrant pairs.
    //**************************************************************************

    // Calculate leakages for the balance table.

    problem_boundary_leakage( data, lkgs_l, fh_i, fv_i, fz_i );

    // Cleanup

    delete [] phij;
    delete [] phik;

    return;
}

void Sweep3d::octant_ordering( Input_edit &data )
{
    // Set up order in which octants will be swept. The order allows parallel
    // sweeps to be ordered across (j,k) for stacks of i-values. This allows
    // us to match the old vectorized sweepers' logic and order which was
    // optimized for the fact that reflective boundaries usually occured at the
    // left(i), bottom(j), and front(k) faces of most problems. Starting sweeps
    // toward reflective boundaries, followed by sweeps away, reduces the amount
    // of angular flux information that must be stored and gives you explicit
    // boundaries rather than implicit ones. jlow, jhigh, klow, khigh, specify
    // the spatial sweep order within an i-plane (a plane normal to x-axis).
    // During loop over iop=1,maxop where maxop=2*mm*it, mu alternates between
    // -1 and 1 rather than doing all mm angles for mu=-1 then all mm angles for
    // mu=1. This allows for reuse of i-direction's angular flux array (phii)
    // at a reflective boundary, which may be completely contained in cache when
    // mu changes sign and the calculation proceeds back along the same angle m.
    // More important than the cache affect, is the avoidance of having to
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

    i_cm   = data.it() + 1;
    ishift = -1;
    ioff   = 0;
    mz     = 0;
    ilim   = 1;

    if ( ih == 0 )
    {
        jshift = -1;
        jlow   = data.jt()-1;
        jhigh  = 0;
    }
    else
    {
        jshift = +1;
        jlow   = 0;
        jhigh  = data.jt()-1;
    }

    if ( k1 == 0 )
    {
        kshift = -1;
        klow   = data.kt()-1;
        khigh  = 0;
    }
    else
    {
        kshift = +1;
        klow   = 0;
        khigh  = data.kt()-1;
    }
}

void Sweep3d::build_angular_source( Input_edit &data, Sn_constants &sn,
                                    Array4D &src_mom, Array3D &phi      )
{
    // Build angular source in phi. phi must first be re-initialized to zero.

    int nk1;  // used to set first index of phi when building
    int nk2;  //   the angular source and the flux moments
    int i;    // loop variable for cells in the x-direction
    int j;    // loop variable for cells in the y-direction
    int k;    // loop variable for cells in the z-direction
    int m;    // loop variable for quadrature points (angles) per quadrant
    int n;    // loop variable for flux and source angular moments

    REAL c1;  // used to hold a value calculated outside of a
    REAL c2;  //   loop, for use inside the loop (efficiency)

    phi.Array3D_reinit( 0.0 );

    for ( n=0 ; n < data.nm() ; n++ )
        for ( m=0 ; m < data.mm() ; m++ )
        {
            nk1 = m * data.it();
            nk2 = nk1 + data.mm() * data.it();
            c1 = sn.p(n,m,iq);
            c2 = sn.p(n,m,iqp);

            for ( i=0 ; i < data.it() ; i++ )
                for ( k=0 ; k < data.kt() ; k++ )
                    for ( j=0 ; j < data.jt() ; j++ )
                    {
                        phi(j,k,i+nk1) += src_mom(j,k,i,n) * c1;
                        phi(j,k,i+nk2) += src_mom(j,k,i,n) * c2;
                    }

        }
}

void Sweep3d::edge_and_boundary_set( Input_edit &data, Array2D &phii,
                                     REAL *phij,       REAL *phik,
                                     Array2D &bsavv,   Array3D &bsavz )
{
    // Initialize edge fluxes, impose boundary conditions. Assumed here, that
    // right, top, and back boundaries are vacuum boundaries. If any of left,
    // bottom, or front boundaries are also vacuum, then their edge flux
    // arrays are set to zero at the time of sweep reversal. If the bottom or
    // front boundary are reflective than phij or phik must be loaded from
    // bsavv or basavz at the time of sweep reversal.

    int j;  // loop variable for cells in the y-direction
    int k;  // loop variable for cells in the z-direction

    if ( ishift == -1 && isw == data.it()-1 )
        phii.Array2D_reinit( 0.0 );
    else if ( ishift == +1 && isw == 0 && data.ibl() == 0 )
        phii.Array2D_reinit( 0.0 );

    if ( ih == 0 )                            // corners 1,2,5,6
        for ( k=0 ; k < data.kt() ; k++ )
            phij[k] = 0.0;
    else                                      // corners 3,4,7,8
    {
        if ( data.ibb() == 0 )
            for ( k=0 ; k < data.kt() ; k++ )
                phij[k] = 0.0;
        else if ( data.ibb() == 1 )
            for ( k=0 ; k < data.kt() ; k++ )
                phij[k] = bsavv(k,iop);
    }

    if ( k1 == 0 )                            // corners 1,2,3,4
        for ( j=0 ; j < data.jt() ; j++ )
            phik[j] = 0.0;
    else                                      // corners 5,6,7,8
    {
        if ( data.ibfr() == 0 )
            for ( j=0 ; j < data.jt() ; j++ )
                phik[j] = 0.0;
        else if ( data.ibfr() == 1 )
            for ( j=0 ; j < data.jt() ; j++ )
                phik[j] = bsavz(j,iop,ih);
    }
}

void Sweep3d::balance_eqn_no_fixup( Input_edit &data, Sn_constants &sn,
                                    Pre_calcs &pre,   Array3D &phi,
                                    Array2D &phii,    REAL *phij,
                                    REAL *phik                          )
{
    // Balance equation with diamond difference, no flux fixup. phii, phij, phik
    // are only needed from one i-plane to the next and thus can be continuously
    // overwritten. phij and phik are saved/restored to/from bsavv and bsavz at
    // boundaries. phii is immediately reused via sweep reversal, thus does not
    // require a storage array.

    int j;    // loop variable for cells in the y-direction
    int k;    // loop variable for cells in the z-direction

    REAL ql;  // a temp variable to hold calculated results

    for ( k=klow ; k != (khigh+kshift) ; k += kshift )
        for ( j=jlow ; j != (jhigh+jshift) ; j += jshift )
        {
            ql = phi(j,k,iz) +
                 phii(j,k) * pre.muh(isw,mz) +
                 phij[k]   * pre.etah(j,mz) +
                 phik[j]   * pre.tsih(k,mz);

            phi(j,k,iz) = ql * pre.dlinv(j,k,isw,mz);

            phii(j,k) = 2.0 * phi(j,k,iz) - phii(j,k);
            phij[k]   = 2.0 * phi(j,k,iz) - phij[k];
            phik[j]   = 2.0 * phi(j,k,iz) - phik[j];
        }
}

void Sweep3d::balance_eqn_with_fixup( Input_edit &data, Sn_constants &sn,
                                      Pre_calcs &pre,   Cross_section &xsec,
                                      Array3D &phi,     Array2D &phii,
                                      REAL *phij,       REAL *phik           )
{
    // Balance equation w/ diamond difference and set to zero flux fixup.

    int j;     // loop variable for cells in the y-direction
    int k;     // loop variable for cells in the z-direction

    REAL ql;      // a temp variable to hold calculated results
    REAL dl_loc;  // a temp variable to hold calculated results
    REAL sih;     // temp holder for phii during set-to-zero fixup
    REAL siv;     // temp holder for phij during set-to-zero fixup
    REAL sif;     // temp holder for phik during set-to-zero fixup

    for ( k=klow ; k != (khigh+kshift) ; k += kshift )
        for ( j=jlow ; j != (jhigh+jshift) ; j += jshift )
        {
            ql = phi(j,k,iz) +
                 phii(j,k) * pre.muh(isw,mz) +
                 phij[k]   * pre.etah(j,mz) +
                 phik[j]   * pre.tsih(k,mz);

            dl_loc = pre.dl(j,k,isw,mz);

            phi(j,k,iz) = ql * pre.dlinv(j,k,isw,mz);

            sih = 2.0 * phi(j,k,iz) - phii(j,k);
            siv = 2.0 * phi(j,k,iz) - phij[k];
            sif = 2.0 * phi(j,k,iz) - phik[j];

            // Test for negative fluxes, if found, adjust temporary values
            // ql and dl, to account for setting the negative flux to zero.
            // After setting a flux to zero and adjusting other fluxes to
            // account for it, go back to make sure the change did not cause
            // other fluxes to become negative.

          neg_flux:

            if ( sih < 0.0 )
            {
                ql -= 0.5 * phii(j,k) * pre.muh(isw,mz);
                dl_loc -= pre.muh(isw,mz);

                phi(j,k,iz) = ql / dl_loc;

                sih = 0.0;
                if ( siv != 0.0 ) siv = 2.0 * phi(j,k,iz) - phij[k];
                if ( sif != 0.0 ) sif = 2.0 * phi(j,k,iz) - phik[j];
            }

            if ( siv < 0.0 )
            {
                ql -= 0.5 * phij[k] * pre.etah(j,mz);
                dl_loc -= pre.etah(j,mz);

                phi(j,k,iz) = ql / dl_loc;

                siv = 0.0;
                if ( sif != 0.0 ) sif = 2.0 * phi(j,k,iz) - phik[j];
                if ( sih != 0.0 ) sih = 2.0 * phi(j,k,iz) - phii(j,k);

                goto neg_flux;
            }

            if ( sif < 0.0 )
            {
                ql -= 0.5 * phik[j] * pre.tsih(k,mz);
                dl_loc -= pre.tsih(k,mz);

                phi(j,k,iz) = ql / dl_loc;

                sif = 0.0;
                if ( sih != 0.0 ) sih = 2.0 * phi(j,k,iz) - phii(j,k);
                if ( siv != 0.0 ) siv = 2.0 * phi(j,k,iz) - phij[k];

                goto neg_flux;
            }

            phii(j,k) = sih;
            phij[k]   = siv;
            phik[j]   = sif;

        }
}

void Sweep3d::save_boundary( Input_edit &data, Array2D &bsavv, Array3D &bsavz,
                             REAL *phij,       REAL *phik                      )
{
    // Save j and k edge fluxes at reflective boundaries. The i-edge fluxes do
    // not need to be saved, they are about to be used in the sweep reversal.

    int j;  // loop variable for cells in the y-direction
    int k;  // loop variable for cells in the z-direction

    if ( data.ibb()  == 1 && ih == 0 )
        for ( k=klow ; k != (khigh+kshift) ; k += kshift )
            bsavv(k,iop) = phij[k];

    if ( data.ibfr() == 1 && k1 == 0 )
        for ( j=jlow ; j != (jhigh+jshift) ; j += jshift )
            bsavz(j,iop,ih) = phik[j];
}

void Sweep3d::cell_boundary_leakage( Input_edit &data, Sn_constants &sn,
                                     Array3D &fh_i,    Array3D &fv_i,
                                     Array3D &fz_i,    Array2D &phii,
                                     REAL *phij,       REAL *phik        )
{
    // Compute leakages.

    int j;    // loop variable for cells in the y-direction
    int k;    // loop variable for cells in the z-direction

    REAL c1;  // used to hold a value calculated outside of a
              //   loop, for use inside the loop (efficiency)

    if ( isw+ioff == 0 && data.ibl() == 0 )
    {
        c1 = ishift * sn.wmu(mz);

        for ( k=0 ; k < data.kt() ; k++ )
            for ( j=0 ; j < data.jt() ; j++ )
                fh_i(j,k,0) += c1 * phii(j,k);
    }

    if ( isw+ioff == data.it() ) 
    {
        c1 = ishift * sn.wmu(mz);

        for ( k=0 ; k < data.kt() ; k++ )
            for ( j=0 ; j < data.jt() ; j++ )
                fh_i(j,k,1) += c1 * phii(j,k);
    }

    if ( ih == 1 || data.ibb() == 0 )
    {
        c1 = jshift * sn.weta(mz);

        for ( k=0 ; k < data.kt() ; k++ )
            fv_i(k,isw,ih) += c1 * phij[k];
    }

    if ( k1 == 1 || data.ibfr() == 0 )
    {
        c1 = kshift * sn.wtsi(mz);

        for ( j=0 ; j < data.jt() ; j++ )
            fz_i(j,isw,k1) += c1 * phik[j];
    }
}

void Sweep3d::sweep_reversal( Input_edit &data )
{
    // Sweep reversal for mu direction cosine.

      if ( i_cm == ilim )
      {
          if ( ishift == -1 )
          {
              i_cm   = 0;
              ishift = +1;
              ioff   = +1;
              ilim   = data.it();
          }
          else if ( mz < data.mm()-1 )
          {
              i_cm   = data.it() + 1;
              ishift = -1;
              ioff   = 0;
              ilim   = 1;
              mz    += 1;
          }
      }
}

void Sweep3d::flux_moments( Input_edit &data, Sn_constants &sn, Array4D &flux,
                            Array3D &phi                                       )
{
    // Compute the flux moments.

    int nk1;  // used to set first index of phi when building
    int nk2;  //   the angular source and the flux moments
    int i;    // loop variable for cells in the x-direction
    int j;    // loop variable for cells in the y-direction
    int k;    // loop variable for cells in the z-direction
    int m;    // loop variable for quadrature points (angles) per quadrant
    int n;    // loop variable for flux and source angular moments

    REAL c1;  // used to hold a value calculated outside of a
    REAL c2;  //   loop, for use inside the loop (efficiency)

    for ( m=0 ; m < data.mm() ; m++ )
    {
        nk1 = m * data.it();
        nk2 = nk1 + data.itmm();

        for ( i=0 ; i < data.it() ; i++ )
            for ( k=0 ; k < data.kt() ; k++ )
                for ( j=0 ; j < data.jt() ; j++ )
                    flux(j,k,i,0) += sn.w(m)*(phi(j,k,i+nk1) + phi(j,k,i+nk2));
    }

    for ( n=1 ; n < data.nm() ; n++ )
        for ( m=0 ; m < data.mm() ; m++ )
        {
            nk1 = m * data.it();
            nk2 = nk1 + data.itmm();
            c1  = sn.w(m) * sn.p(n,m,iq);
            c2  = sn.w(m) * sn.p(n,m,iqp);

            for ( i=0 ; i < data.it() ; i++ )
                for ( k=0 ; k < data.kt() ; k++ )
                    for ( j=0 ; j < data.jt() ; j++ )
                        flux(j,k,i,n) += c1*phi(j,k,i+nk1) + c2*phi(j,k,i+nk2);
        }
}

void Sweep3d::problem_boundary_leakage( Input_edit &data, REAL *lkgs_l,
                                        Array3D &fh_i,    Array3D &fv_i,
                                        Array3D &fz_i                    )
{
    // Calculate leakages for the balance table.

    int i;  // loop variable for cells in the x-direction
    int j;  // loop variable for cells in the y-direction
    int k;  // loop variable for cells in the z-direction

    lkgs_l[0] = 0.0;  // left+right leakage (horizontal or x-axis)
    lkgs_l[1] = 0.0;  // right      leakage
    lkgs_l[2] = 0.0;  // bottom+top leakage (vertical or y-axis)
    lkgs_l[3] = 0.0;  // top        leakage
    lkgs_l[4] = 0.0;  // front+back leakage (in-out or z-axis)
    lkgs_l[5] = 0.0;  // back       leakage

    for ( k=0 ; k < data.kt() ; k++ )
        for ( j=0 ; j < data.jt() ; j++ )
        {
            fh_i(j,k,0) *= 4.0 / (data.hj(j) * data.hk(k));
            fh_i(j,k,1) *= 4.0 / (data.hj(j) * data.hk(k));
            lkgs_l[0] -= fh_i(j,k,0);
            lkgs_l[1] += fh_i(j,k,1);
        }

    for ( i=0 ; i < data.it() ; i++ )
        for ( k=0 ; k < data.kt() ; k++ )
        {
            fv_i(k,i,0) *= 4.0 / (data.hi(i) * data.hk(k));
            fv_i(k,i,1) *= 4.0 / (data.hi(i) * data.hk(k));
            lkgs_l[2] -= fv_i(k,i,0);
            lkgs_l[3] += fv_i(k,i,1);
        } 

    for ( i=0 ; i < data.it() ; i++ )
        for ( j=0 ; j < data.jt() ; j++ )
        {
            fz_i(j,i,0) *= 4.0 / (data.hi(i) * data.hj(j));
            fz_i(j,i,1) *= 4.0 / (data.hi(i) * data.hj(j));
            lkgs_l[4] -= fz_i(j,i,0);
            lkgs_l[5] += fz_i(j,i,1);
        }

    lkgs_l[0] += lkgs_l[1];
    lkgs_l[2] += lkgs_l[3];
    lkgs_l[4] += lkgs_l[5];
}

//---------------------------------------------------------------------------//
//                              end of Sweep3d.cc
//---------------------------------------------------------------------------//

