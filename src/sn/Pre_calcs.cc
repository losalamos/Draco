//----------------------------------*-C++-*----------------------------------//
// Pre_calcs.cc
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Pre-calculate all possible factors for use in inner solver loops.
//---------------------------------------------------------------------------//

#include "sn/Pre_calcs.hh"

Pre_calcs::Pre_calcs( Input_edit &data, Cross_section &xsec, Sn_constants &sn ):
             dlinv_p( data.jt(), data.kt(), data.it(), data.mm() ),
             dl_p( data.jt(), data.kt(), data.it(), data.mm() ),
             muh_p( data.it(), data.mm() ), etah_p( data.jt(), data.mm() ),
             tsih_p( data.kt(), data.mm() )
{
    for ( int m=0 ; m < data.mm() ; m++ )
    {
        for ( int i=0 ; i < data.it() ; i++ )
            muh_p(i,m) = sn.mu(m)  * data.hi(i);
        for ( int j=0 ; j < data.jt() ; j++ )
            etah_p(j,m) = sn.eta(m) * data.hj(j);
        for ( int k=0 ; k < data.kt() ; k++ )
            tsih_p(k,m) = sn.tsi(m) * data.hk(k);
    }

    for ( int m=0 ; m < data.mm() ; m++ )
        for ( int i=0 ; i < data.it() ; i++ )
            for ( int k=0 ; k < data.kt() ; k++ )
                for ( int j=0 ; j < data.jt() ; j++ )
                    dlinv_p(j,k,i,m) = 1.0 / ( xsec.ct(j,k,i) + muh_p(i,m) +
                                               etah_p(j,m) + tsih_p(k,m)     );

    for ( int m=0 ; m < data.mm() ; m++ )
        for ( int i=0 ; i < data.it() ; i++ )
            for ( int k=0 ; k < data.kt() ; k++ )
                for ( int j=0 ; j < data.jt() ; j++ )
                    dl_p(j,k,i,m) = xsec.ct(j,k,i) + muh_p(i,m) +
                                    etah_p(j,m) + tsih_p(k,m);
}

//---------------------------------------------------------------------------//
//                              end of Pre_calcs.cc
//---------------------------------------------------------------------------//

