//----------------------------------*-C++-*----------------------------------//
// Test_3T.cc
// Geoffrey M. Furnish
// Wed Nov 19 17:05:08 1997
//---------------------------------------------------------------------------//
// @> Test problem template.
//---------------------------------------------------------------------------//

#include "3T/test/Test_3T.hh"

#include <iostream.h>

template<class Problem>
Test_3T<Problem>::Test_3T()
    : Problem()
{
}


template<class Problem>
void Test_3T<Problem>::run()
{
    cout << "Running the test problem, NOT!" << endl;
}

//---------------------------------------------------------------------------//
//                              end of Test_3T.cc
//---------------------------------------------------------------------------//
