//----------------------------------*-C++-*----------------------------------//
// Inner_iter.hh
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Controller of inner iterations.
//---------------------------------------------------------------------------//

#ifndef __sn_Inner_iter_hh__
#define __sn_Inner_iter_hh__

#include "sn/precision.hh"
#include "sn/Array.hh"
#include "sn/Error.hh"
#include "sn/Input_edit.hh"
#include "sn/Output_edit.hh"
#include "sn/Source.hh"
#include "sn/Sweep3d.hh"

class Inner_iter
{

    public:

        // use the default constructor
        // Inner_iter();

        // use the default destructor
        // ~Inner_iter();

        void do_inner_iter( Input_edit &data );

    private:

};

#endif                          // __sn_Inner_iter_hh__

//---------------------------------------------------------------------------//
//                              end of Inner_iter.hh
//---------------------------------------------------------------------------//

