//----------------------------------*-C++-*----------------------------------//
// Pre_calcs.cc
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Pre-calculate all possible factors for use in inner solver loops.
//---------------------------------------------------------------------------//

#include "sn/Pre_calcs.hh"

Pre_calcs::Pre_calcs( Input_edit &data, Cross_section &xsec, Sn_constants &sn ):
             dlinv_p( data.jt(), data.kt(), data.it(), data.mm() )
{
    for ( int m=0 ; m < data.mm() ; m++ )
        for ( int i=0 ; i < data.it() ; i++ )
            for ( int k=0 ; k < data.kt() ; k++ )
                for ( int j=0 ; j < data.jt() ; j++ )
                    dlinv_p(j,k,i,m) = 1.0 / ( xsec.ct(j,k,i) +
                                               sn.mu(m)  * data.hi(i) +
                                               sn.eta(m) * data.hj(j) +
                                               sn.tsi(m) * data.hk(k)   );
}

//---------------------------------------------------------------------------//
//                              end of Pre_calcs.cc
//---------------------------------------------------------------------------//

