//----------------------------------*-C++-*----------------------------------//
// Test_3T.cc
// Geoffrey M. Furnish
// Wed Nov 19 17:05:08 1997
//---------------------------------------------------------------------------//
// @> Test problem template.
//---------------------------------------------------------------------------//

#include "3T/test/Test_3T.hh"

#include <iostream.h>

template<class MT, class Problem>
Test_3T<MT, Problem>::Test_3T( const SP<MT>& spm_, const Quad_Params& q )
    : Problem()
{
}


template<class MT, class Problem>
void Test_3T<MT, Problem>::run()
{
    cout << "Running the test problem, NOT!" << endl;
}

//---------------------------------------------------------------------------//
//                              end of Test_3T.cc
//---------------------------------------------------------------------------//
