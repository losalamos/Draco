//----------------------------------*-C++-*----------------------------------//
// Cross_section.hh
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Define cross sections for the problem.
//---------------------------------------------------------------------------//

#ifndef __sn_Cross_section_hh__
#define __sn_Cross_section_hh__

#include "sn/precision.hh"
#include "sn/Array.hh"
#include "sn/Input_edit.hh"

class Cross_section
{

    public:

        // the constructor is used to set the cross section values

        Cross_section( Input_edit &data );

        // use the default destructor
        // ~Cross_section();

        REAL ct( int j, int k, int i ) const
        {
            return ct_p(j,k,i);
        }

        REAL sigs( int j, int k, int i, int isct_l ) const
        {
            return sigs_p(j,k,i,isct_l);
        }

    private:

        Array3D ct_p;    // total cross section
        Array4D sigs_p;  // macroscopic differential scattering
                         //   cross section

};

#endif                          // __sn_Cross_section_hh__

//---------------------------------------------------------------------------//
//                              end of Cross_section.hh
//---------------------------------------------------------------------------//

