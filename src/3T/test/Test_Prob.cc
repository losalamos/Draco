//----------------------------------*-C++-*----------------------------------//
// Test_Prob.cc
// Geoffrey M. Furnish
// Wed Nov 19 16:18:54 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/test/Test_Prob.hh"
#include "3T/test/Test_3T.hh"
#include "3T/test/XYZ_Quadratic.hh"

SP<Test_Prob> Test_Prob_allocator( int argc, char *argv[] )
{
    SP<Test_Prob> prob;

// Theoretically we could parse argc, argv to figure out which test problem
// to initiate.  For now, however, we just hardwire one.

    prob = new Test_3T<XYZ_Quadratic>();

    return prob;
}

//---------------------------------------------------------------------------//
//                              end of Test_Prob.cc
//---------------------------------------------------------------------------//
