//----------------------------------*-C++-*----------------------------------//
// Error.cc
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Error convergence.
//---------------------------------------------------------------------------//

#include "sn/Error.hh"

#include <iostream.h>

Error::Error( Input_edit &data ) : pflux(data.jt(),data.kt(),data.it())
{}

void Error::error_init( Input_edit &data, Array4D &flux )
{
    // save the flux from the previous iteration

    for ( int i=0 ; i < data.it() ; i++ )
        for ( int k=0 ; k < data.kt() ; k++ )
            for ( int j=0 ; j < data.jt() ; j++ )
                pflux(j,k,i) = flux(j,k,i,0);
}

void Error::error_calc( Input_edit &data, Array4D &flux )
{
    // store the absolute relative difference between the current and previous
    // fluxes in pflux (in other words replace pflux values with relative error)

    for ( int i=0 ; i < data.it() ; i++ )
        for ( int k=0 ; k < data.kt() ; k++ )
            for ( int j=0 ; j < data.jt() ; j++ )
                if ( flux(j,k,i,0) != 0.0 )
                {
                    pflux(j,k,i) = ( flux(j,k,i,0) - pflux(j,k,i) ) /
                                     flux(j,k,i,0);
                    pflux(j,k,i) = ( pflux(j,k,i) > 0.0 ) ? pflux(j,k,i) :
                                                          - pflux(j,k,i) ;
                }
                else
                    pflux(j,k,i) = 0.0;

    // find the maximum relative error

    err = 0.0;

    for ( int i=0 ; i < data.it() ; i++ )
        for ( int k=0 ; k < data.kt() ; k++ )
            for ( int j=0 ; j < data.jt() ; j++ )
                if ( err < pflux(j,k,i) )
                     err = pflux(j,k,i);
}

void Error::error_print( int its )
{
    cout << "its = " << its << "  err = " << err << endl;
}

bool Error::error_test( Input_edit &data, int its )
{
    bool converged = 0;

    if ( ( err <= data.epsi() && data.epsi() >= 0.0 ) ||
         ( its >= int(-data.epsi()+.99) && data.epsi() < 0.0 ) )
        converged = 1;

    return converged;
}

//---------------------------------------------------------------------------//
//                              end of Error.cc
//---------------------------------------------------------------------------//

