//----------------------------------*-C++-*----------------------------------//
// Pre_calcs.hh
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Pre-calculate all possible factors for use in inner solver loops.
//---------------------------------------------------------------------------//

#ifndef __sn_Pre_calcs_hh__
#define __sn_Pre_calcs_hh__

#include "sn/precision.hh"
#include "sn/Array.hh"
#include "sn/Cross_section.hh"
#include "sn/Input_edit.hh"
#include "sn/Sn_constants.hh"

class Pre_calcs
{

    public:

        // the constructor calculates all factors

        Pre_calcs( Input_edit &data, Cross_section &xsec, Sn_constants &sn );

        // use the default destructor
        // ~Pre_calcs();

        REAL dlinv( int j, int k, int i, int m ) const
        {
            return dlinv_p(j,k,i,m);
        }

    private:

        Array4D dlinv_p;  // a factor in the balance equation

};

#endif                          // __sn_Pre_calcs_hh__

//---------------------------------------------------------------------------//
//                              end of Pre_calcs.hh
//---------------------------------------------------------------------------//

