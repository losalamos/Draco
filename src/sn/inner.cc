//----------------------------------*-C++-*----------------------------------//
// inner.cc
// Scott Turner
// 19 February 1998
//---------------------------------------------------------------------------//
// @> Controller of inner iterations.
//---------------------------------------------------------------------------//

//
// Inner sets up cross sections and geometry for the inner iteration sweeps.
// It also tests for convergence, does timing and prints results
//

#include <math.h>
#include <iostream.h>
#include <stdlib.h>

#include "sn/test/snpp.hh"
#include "sn/test/protos.hh"

void inner( int it,     int jt,    int kt,     int mm,    int nm,
            int isct,   int isctp, int ibl,    int ibb,   int ibfr,
            int iprint, int ifxg,  REAL dx,    REAL dy,   REAL dz,
            REAL epsi                                              )
{

#include "sn/test/array.hh"

  int its;      // convergence iteration counter
  int s1;       // used for setting up the spherical harmonics
  int s2;       // used for setting up the spherical harmonics
  int s3;       // used for setting up the spherical harmonics
  int i;        // loop variable for cells in the x-direction
  int j;        // loop variable for cells in the y-direction
  int k;        // loop variable for cells in the z-direction
  int m;        // loop variable for quadrature points (angles) per quadrant
  int n;        // loop variable for flux and source angular moments
  int iq;       // loop variable for the number of quadrants, may be calculated,
                //   rather than incremented in some cases

  REAL cpu0;    // starting time for determining total cpu       time
  REAL cpu1;    // ending   time for determining total cpu       time
  REAL wall0;   // starting time for determining total wallclock time
  REAL wall1;   // ending   time for determining total wallclock time
  REAL grindw;  // wallclock grind time
  REAL grindc;  // cpu       grind time
  REAL fsum;    // scalar flux sum over the entire mesh
  REAL err;     // max change in the flux from one iteration to the next

  Array3D p(nm,mm,8);           // spherical harmonics array
  Array3D ct(jt,kt,it);         // the total cross section in each mesh cell
  Array4D src(jt,kt,it,nm);     // the source in each cell for each moment
  Array4D flux(jt,kt,it,nm);    // the scalar flux in each cell for each moment
  Array3D pflux(jt,kt,it);      // both the previous flux and the relative
                                //   change in the flux after each iteration
  Array3D srcx(jt,kt,it);       // fixed source
  Array4D sigs(jt,kt,it,isctp); // macroscopic differential scattering
                                //   cross section

// allocate memory for local arrays

  REAL *const lkgs_l = new REAL [6];  // the leakages from each face of the cube
  REAL *const w      = new REAL [mm]; // sn constants
  REAL *const mu     = new REAL [mm]; // sn constants
  REAL *const eta    = new REAL [mm]; // sn constants
  REAL *const tsi    = new REAL [mm]; // sn constants
  REAL *const wmu    = new REAL [mm]; // sn constants
  REAL *const weta   = new REAL [mm]; // sn constants
  REAL *const wtsi   = new REAL [mm]; // sn constants

// set precision for output of floating point numbers

  cout.precision(13);

// initialize cross sections and source

  for ( i=0 ; i < it ; i++ ) {
  for ( k=0 ; k < kt ; k++ ) {
  for ( j=0 ; j < jt ; j++ ) {
    srcx(j,k,i)     = 1.0;    // uniform fixed source
    ct(j,k,i)       = 1.0;    // uniform total cross section
    sigs(j,k,i,0)   = 0.5;    // scattering ratio (c) = 0.5
    if (isct == 1)
      sigs(j,k,i,1) = 0.6;    // linear anisotropic (2l+1)mubar, l=1
  } } }

// define the Sn constants

  if (mm == 6)                      // S6
  {    
    mu[0]  = 0.23009194;
    eta[0] = 0.94557676;
    w[0]   = 0.16944656 / 8.0;
    
    mu[1]  = 0.68813432;
    eta[1] = mu[1];
    w[1]   = 0.16388677 / 8.0;
  
    mu[2]  = mu[0];
    eta[2] = mu[1];
    w[2]   = w[1];
  
    mu[3]  = eta[0]; 
    eta[3] = mu[0];
    w[3]   = w[0];
  
    mu[4]  = mu[1];
    eta[4] = mu[0];
    w[4]   = w[1];
  
    mu[5]  = mu[0];
    eta[5] = mu[0];
    w[5]   = w[0];
  }
  else if (mm == 3)                 // S4
  {
    mu[0]  = 0.30163878;
    eta[0] = 0.90444905;
    w[0]   = 1.0 / 3.0 / 8.0;
  
    mu[1]  = eta[0];
    eta[1] = mu[0];
    w[1]   = w[0];
  
    mu[2]  = mu[0];
    eta[2] = mu[0];
    w[2]   = w[0];
  }

  for ( m = 0 ; m < mm ; m++ )
  {
    tsi[m]  = sqrt(1.-mu[m]*mu[m]-eta[m]*eta[m]);
    wmu[m]  = w[m] * mu[m];
    weta[m] = w[m] * eta[m];
    wtsi[m] = w[m] * tsi[m];
  }

// spherical harmonics

  for ( m  = 0 ; m  < mm ;  m++ )  {
  for ( iq = 0 ; iq < 8  ; iq++ )
    p(0,m,iq) = 1.0;
  }

  if (isct > 0)                          {
      iq = -1;
    for ( s1 = -1 ; s1 < 2 ; s1 += 2 )   { 
    for ( s2 = -1 ; s2 < 2 ; s2 += 2 )   { 
    for ( s3 = -1 ; s3 < 2 ; s3 += 2 )   { 
      ++iq;
      for ( m  = 0 ; m  < mm ;  m++ )    {
        p(1,m,iq) = s3 * mu[m];
        p(2,m,iq) = s2 * eta[m];
        p(3,m,iq) = s1 * tsi[m];
  } } } } }

// initialize all flux moments and the iteration counter to zero
      
  flux.Array4D_reinit(0.0);

  its = 0;

// begin cpu and wallclock timing

  timer( cpu0, wall0 );

// begin iterations

  while (1)
  {

    its += 1;

// build the source, first add the fixed source to the isotropic scatter source

    for ( i=0 ; i < it ; i++ ) {
    for ( k=0 ; k < kt ; k++ ) {
    for ( j=0 ; j < jt ; j++ ) {
      src(j,k,i,0) = srcx(j,k,i) + sigs(j,k,i,0) * flux(j,k,i,0);
    } } }

// now, compute the (linearly) anisotropic scattering components (also referred
// to as higher source moments)

    for ( n=1 ; n < nm ; n++ ) {
    for ( i=0 ; i < it ; i++ ) {
    for ( k=0 ; k < kt ; k++ ) {
    for ( j=0 ; j < jt ; j++ ) {
      src(j,k,i,n) = sigs(j,k,i,1) * flux(j,k,i,n);
    } } } }

// save the previous scalar flux for use in the convergence test, and
// re-initialize the current flux moments to zero

    for ( i=0 ; i < it ; i++ ) {
    for ( k=0 ; k < kt ; k++ ) {
    for ( j=0 ; j < jt ; j++ ) {
      pflux(j,k,i) = flux(j,k,i,0);
    } } }

    flux.Array4D_reinit(0.0);

// perform an ordered sweep through the mesh.

   sweep3d ( it, jt, kt, mm, nm,
             ibl,   ibb,     ibfr,
             ifxg,  dx,      dy,
             dz,    lkgs_l, w,
             mu,    eta,    tsi,
             wmu,   weta,   wtsi,
             p,     ct,     src,
             flux                  );

// store the absolute relative difference between the current and previous
// fluxes in pflux (in other words replace pflux values with relative error).

    for ( i=0 ; i < it ; i++ ) {
    for ( k=0 ; k < kt ; k++ ) {
    for ( j=0 ; j < jt ; j++ ) {
      if (flux(j,k,i,0) != 0.0)
      {
        pflux(j,k,i) = (flux(j,k,i,0) - pflux(j,k,i)) /
                        flux(j,k,i,0);
        pflux(j,k,i) = ( pflux(j,k,i) > 0.0 ) ? pflux(j,k,i) : -pflux(j,k,i);
      }
      else
        pflux(j,k,i) = 0.0;
    } } }

// find the maximum relative error

    err = 0.0;

    for ( i=0 ; i < it ; i++ ) {
    for ( k=0 ; k < kt ; k++ ) {
    for ( j=0 ; j < jt ; j++ ) {
      if ( err < pflux(j,k,i) ) err = pflux(j,k,i);
    } } }
      
// print iteration/error information and test for convergence

    cout << "its = " << its << "  err = " << err << endl;

    if ( (err <= epsi && epsi >= 0.0) ||
         (its >= int(-epsi+.99) && epsi < 0.0) )

      break;

  }

// stop cpu and wallclock timing

  timer( cpu1, wall1 );

// sum the scalar fluxes over the entire mesh

  fsum = 0.0;

  for ( i=0 ; i < it ; i++ ) {
  for ( k=0 ; k < kt ; k++ ) {
  for ( j=0 ; j < jt ; j++ ) {
    fsum += flux(j,k,i,0);
  } } }

// print out the total scalar flux and the leakages from each face

  cout << endl << "scalar flux sum = " << fsum << endl << endl;
  cout << "left + right leakage = " << lkgs_l[0] << endl;
  cout << "       right leakage = " << lkgs_l[1] << endl;
  cout << "bottom + top leakage = " << lkgs_l[2] << endl;
  cout << "         top leakage = " << lkgs_l[3] << endl;
  cout << "front + back leakage = " << lkgs_l[4] << endl;
  cout << "        back leakage = " << lkgs_l[5] << endl << endl;

// calculate grind times and print them out with the cpu and wallclock times
      
  grindw = 1.0e9 * (wall1-wall0) / ( REAL(it)*REAL(jt)*REAL(kt)*
                                     REAL(mm)*REAL(8)*REAL(its)  );
  grindc = 1.0e9 * (cpu1 - cpu0) / ( REAL(it)*REAL(jt)*REAL(kt)*
                                     REAL(mm)*REAL(8)*REAL(its)  );
  cout << " CPU       time for Sn: " << cpu1 - cpu0 << endl;
  cout << " Wallclock time for Sn: " << wall1-wall0 << endl;
  cout << " CPU        grind time: " << grindc      << endl;
  cout << " Wallclock  grind time: " << grindw      << endl;

  if (iprint >= 1)                             {
    for ( n=0 ; n < nm ; n++ )                 {
    for ( k=0 ; k < kt ; k++ )                 {
      cout << "component number " << n << endl;
      for ( j=jt-1 ; j > -1 ; j-- )            {
      for ( i=0 ; i < it ; i++ )
        cout << " " << flux(j,k,i,n);
      cout << endl;
  } } } }

  delete [] lkgs_l;
  delete [] w     ;
  delete [] mu    ;
  delete [] eta   ;
  delete [] tsi   ;
  delete [] wmu   ;
  delete [] weta  ;
  delete [] wtsi  ;

  return;
}

//---------------------------------------------------------------------------//
//                              end of inner.cc
//---------------------------------------------------------------------------//

