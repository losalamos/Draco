//----------------------------------*-C++-*----------------------------------//
// test_utils.cc
// John McGhee
// Thu Aug 27 13:33:43 1998
//---------------------------------------------------------------------------//
// @> Utility function used in the time-step tester.
//---------------------------------------------------------------------------//

#include "test_utils.hh"

#include <cmath>

namespace rtt_linalg
{
 bool compare_reals(const double x1, const double x2, const int ndigits)
 {
     // Determines if two reals are the same to the specified number of
     // significant decimal digits.
    
     double zz;
     if (x2 != 0.)
     {
	 zz = std::log10 ( std::abs( (x1 - x2)/x2 ) );
     }
     else
     {
	 zz = std::log10 ( std::abs( x1 ) );
     }

     return zz < -float(ndigits);

 }

 std::string testMsg(bool testPassed)
 {
     if (testPassed)
	 return "test: passed";
     else
	 return "test: failed";
 }
 
}
//---------------------------------------------------------------------------//
//                              end of test_utils.cc
//---------------------------------------------------------------------------//
