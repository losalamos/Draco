//----------------------------------*-C++-*----------------------------------//
// test_utils.hh
// John McGhee
// Thu Aug 27 13:33:42 1998
//---------------------------------------------------------------------------//
// @> Utility functions used in time-step tester
//---------------------------------------------------------------------------//

#ifndef __timestep_test_test_utils_hh__
#define __timestep_test_test_utils_hh__

#include <string>

namespace rtt_pcgWrap
{
 bool compare_reals(const double x1, const double x2, const int ndigits);
 std::string testMsg(bool testPassed);
}

#endif                          // __timestep_test_test_utils_hh__

//---------------------------------------------------------------------------//
//                              end of test_utils.hh
//---------------------------------------------------------------------------//
