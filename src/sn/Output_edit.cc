//----------------------------------*-C++-*----------------------------------//
// Output_edit.cc
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Define the source.
//---------------------------------------------------------------------------//

#include "sn/Output_edit.hh"

#include <iostream.h>

void Output_edit::print_flux_sum( Input_edit &data, Mat4<REAL> &flux )
{
    // sum the scalar flux over the entire mesh, then print it out

    flux_sum = 0.0;

    for ( int i=0 ; i < data.it() ; i++ )
        for ( int k=0 ; k < data.kt() ; k++ )
            for ( int j=0 ; j < data.jt() ; j++ )
                flux_sum += flux(j,k,i,0);

    cout << endl << "scalar flux sum = " << flux_sum << endl << endl;
}

void Output_edit::print_leakages( Mat1<REAL> &lkgs_l )
{
    cout << "left + right leakage = " << lkgs_l(0) << endl;
    cout << "       right leakage = " << lkgs_l(1) << endl;
    cout << "bottom + top leakage = " << lkgs_l(2) << endl;
    cout << "         top leakage = " << lkgs_l(3) << endl;
    cout << "front + back leakage = " << lkgs_l(4) << endl;
    cout << "        back leakage = " << lkgs_l(5) << endl << endl;
}

void Output_edit::print_timing( REAL cpu0, REAL cpu1, REAL wall0,
                                REAL wall1, Input_edit &data, int its )
{
    // calculate grind times and print them out with the cpu and wallclock times

    REAL grindw;  // wallclock grind time
    REAL grindc;  // cpu       grind time

    grindw = 1.0e9 * ( wall1 - wall0 ) /
             ( REAL(data.it()) * REAL(data.jt()) * REAL(data.kt()) *
               REAL(data.mm()) * REAL(8) * REAL(its)                 );
    grindc = 1.0e9 * ( cpu1 - cpu0 ) /
             ( REAL(data.it()) * REAL(data.jt()) * REAL(data.kt()) *
               REAL(data.mm()) * REAL(8) * REAL(its)                 );

    cout << " CPU       time for Sn: " << cpu1 - cpu0   << endl;
    cout << " Wallclock time for Sn: " << wall1 - wall0 << endl;
    cout << " CPU        grind time: " << grindc        << endl;
    cout << " Wallclock  grind time: " << grindw        << endl;
}

void Output_edit::print_flux_moments( Input_edit &data, Mat4<REAL> &flux )
{
    if ( data.iprint() >= 1 )
    {
        for ( int n=0 ; n < data.nm() ; n++ )
        for ( int k=0 ; k < data.kt() ; k++ )
        {
            cout << "component number " << n << endl;

            for ( int j=data.jt()-1 ; j > -1 ; j-- )
            for ( int i=0 ; i < data.it() ; i++ )
                cout << " " << flux(j,k,i,n);

            cout << endl;
        }
    }
}

//---------------------------------------------------------------------------//
//                              end of Output_edit.cc
//---------------------------------------------------------------------------//

