//----------------------------------*-C++-*----------------------------------//
// Source.hh
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Define the source.
//---------------------------------------------------------------------------//

#ifndef __sn_Source_hh__
#define __sn_Source_hh__

#include "sn/precision.hh"
#include "sn/Array.hh"
#include "sn/Cross_section.hh"
#include "sn/Input_edit.hh"

class Source
{

    public:

        // the constructor is used to set the fixed source

        Source( Input_edit &data );

        // use the default destructor
        // ~Source();

        void build_source( Input_edit &data, Cross_section &xsec,
                           Array4D &src_mom, Array4D &flux        );

    private:

        Array3D fixed_src;

};

#endif                          // __sn_Source_hh__

//---------------------------------------------------------------------------//
//                              end of Source.hh
//---------------------------------------------------------------------------//

