//----------------------------------*-C++-*----------------------------------//
// Cross_section.cc
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Define cross sections for the problem.
//---------------------------------------------------------------------------//

#include "sn/Cross_section.hh"

Cross_section::Cross_section( Input_edit &data ) :
    ct_p(data.jt(),data.kt(),data.it()),
    sigs_p(data.jt(),data.kt(),data.it(),data.isctp())
{
    for ( int i=0 ; i < data.it() ; i++ )
        for ( int k=0 ; k < data.kt() ; k++ )
            for ( int j=0 ; j < data.jt() ; j++ )
            {
                ct_p(j,k,i)     = 1.0;
                sigs_p(j,k,i,0) = 0.5;      // scattering ratio (c) = 0.5

                if ( data.isct() == 1 )
                    sigs_p(j,k,i,1) = 0.6;  // linear anisotropic
                                            // (2l+1)*mubar with l=1, mubar=0.2
            }
}

//---------------------------------------------------------------------------//
//                              end of Cross_section.cc
//---------------------------------------------------------------------------//

