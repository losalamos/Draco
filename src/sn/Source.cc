//----------------------------------*-C++-*----------------------------------//
// Source.cc
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Define the source.
//---------------------------------------------------------------------------//

#include "sn/Source.hh"

Source::Source( Input_edit &data ) : fixed_src(data.jt(),data.kt(),data.it())
{
    for ( int i=0 ; i < data.it() ; i++ )
        for ( int k=0 ; k < data.kt() ; k++ )
            for ( int j=0 ; j < data.jt() ; j++ )
                fixed_src(j,k,i) = 1.0;
}

void Source::build_source( Input_edit &data, Cross_section &xsec,
                           Array4D &src_mom, Array4D &flux        )
{
    // add the fixed source to the isotropic scatter source

    for ( int i=0 ; i < data.it() ; i++ )
        for ( int k=0 ; k < data.kt() ; k++ )
            for ( int j=0 ; j < data.jt() ; j++ )
                src_mom(j,k,i,0) = fixed_src(j,k,i) + xsec.sigs(j,k,i,0) *
                                                      flux(j,k,i,0);

    // now, compute the (linearly) anisotropic scattering components (also
    // referred to as higher source moments)

    for ( int n=1 ; n < data.nm() ; n++ )
        for ( int i=0 ; i < data.it() ; i++ )
            for ( int k=0 ; k < data.kt() ; k++ )
                for ( int j=0 ; j < data.jt() ; j++ )
                    src_mom(j,k,i,n) = xsec.sigs(j,k,i,1) * flux(j,k,i,n);
}

//---------------------------------------------------------------------------//
//                              end of Source.cc
//---------------------------------------------------------------------------//

