//----------------------------------*-C++-*----------------------------------//
// Sn_constants.cc
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Define the Sn constants.
//---------------------------------------------------------------------------//

#include "sn/Sn_constants.hh"

#include <math.h>

Sn_constants::Sn_constants( Input_edit &data ) : p_p(data.nm(),data.mm(),8)
{
    // quadrature

    w_p    = Mat1<REAL>(data.mm());
    mu_p   = Mat1<REAL>(data.mm());
    eta_p  = Mat1<REAL>(data.mm());
    tsi_p  = Mat1<REAL>(data.mm());
    wmu_p  = Mat1<REAL>(data.mm());
    weta_p = Mat1<REAL>(data.mm());
    wtsi_p = Mat1<REAL>(data.mm());

    if ( data.mm() == 6 )           // S6
    {
        mu_p(0)  = 0.23009194;
        eta_p(0) = 0.94557676;
        w_p(0)   = 0.16944656 / 8.0;

        mu_p(1)  = 0.68813432;
        eta_p(1) = mu_p(1);
        w_p(1)   = 0.16388677 / 8.0;

        mu_p(2)  = mu_p(0);
        eta_p(2) = mu_p(1);
        w_p(2)   = w_p(1);
  
        mu_p(3)  = eta_p(0);
        eta_p(3) = mu_p(0);
        w_p(3)   = w_p(0);
  
        mu_p(4)  = mu_p(1);
        eta_p(4) = mu_p(0);
        w_p(4)   = w_p(1);
  
        mu_p(5)  = mu_p(0);
        eta_p(5) = mu_p(0);
        w_p(5)   = w_p(0);
    }
    else if ( data.mm() == 3 )      // S4
    {
        mu_p(0)  = 0.30163878;
        eta_p(0) = 0.90444905;
        w_p(0)   = 1.0 / 3.0 / 8.0;

        mu_p(1)  = eta_p(0);
        eta_p(1) = mu_p(0);
        w_p(1)   = w_p(0);

        mu_p(2)  = mu_p(0);
        eta_p(2) = mu_p(0);
        w_p(2)   = w_p(0);
    }

    for ( int m = 0 ; m < data.mm() ; m++ )
    {
        tsi_p(m)  = sqrt( 1.0 - mu_p(m) * mu_p(m) - eta_p(m) * eta_p(m) );
        wmu_p(m)  = w_p(m) * mu_p(m);
        weta_p(m) = w_p(m) * eta_p(m);
        wtsi_p(m) = w_p(m) * tsi_p(m);
    }

    // spherical harmonics

    int m;
    int iq;

    for ( m = 0 ; m < data.mm() ; m++ )
        for ( iq = 0 ; iq < 8 ; iq++ )
            p_p(0,m,iq) = 1.0;

    if ( data.isct() > 0 )
    {
        iq = -1;

        for ( int s1 = -1 ; s1 < 2 ; s1 += 2 )
            for ( int s2 = -1 ; s2 < 2 ; s2 += 2 )
                for ( int s3 = -1 ; s3 < 2 ; s3 += 2 )
                {
                    ++iq;

                    for ( int m = 0 ; m < data.mm() ;  m++ )
                    {
                        p_p(1,m,iq) = s3 * mu_p(m);
                        p_p(2,m,iq) = s2 * eta_p(m);
                        p_p(3,m,iq) = s1 * tsi_p(m);
                    }
                }
    }
}

//---------------------------------------------------------------------------//
//                              end of Sn_constants.cc
//---------------------------------------------------------------------------//

