//----------------------------------*-C++-*----------------------------------//
// Test_3T.hh
// Geoffrey M. Furnish
// Wed Nov 19 17:05:08 1997
//---------------------------------------------------------------------------//
// @> Test problem template.
//---------------------------------------------------------------------------//

#ifndef __3T_test_Test_3T_hh__
#define __3T_test_Test_3T_hh__

#include "Test_Prob.hh"

//===========================================================================//
// class Test_3T - Template class for accomodating various formulations

// This class implements the Test_Prob abstraction, but is not itself a
// "complete" test.  Rather, it is a template for a test, and the specific
// test problem is provided as a template parameter.
//===========================================================================//

template<class Problem>
class Test_3T : public Test_Prob,
		private Problem
{

  public:
    Test_3T();
//     Test_3T( const Test_3T& );
//     ~Test_3T();
//     Test_3T& operator=( const Test_3T& );

    void run();
};

#endif                          // __3T_test_Test_3T_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/test/Test_3T.hh
//---------------------------------------------------------------------------//
