//----------------------------------*-C++-*----------------------------------//
// Test_Prob.hh
// Geoffrey M. Furnish
// Wed Nov 19 16:18:54 1997
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_test_Test_Prob_hh__
#define __3T_test_Test_Prob_hh__

#include "ds++/SP.hh"

//===========================================================================//
// class Test_Prob - Abstract base class for test problems

// This class represents the abstraction of a test problem for the 3T
// package.  To make a test problem, derive from this, and implement the pure
// virtual methods.
//===========================================================================//

class Test_Prob {

  public:
    virtual ~Test_Prob() {}
    virtual void run() =0;
};

dsxx::SP<Test_Prob> Test_Prob_allocator( int argc, char *argv[] );

#endif                          // __3T/test_Test_Prob_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/test/Test_Prob.hh
//---------------------------------------------------------------------------//
