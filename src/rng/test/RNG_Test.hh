//----------------------------------*-C++-*----------------------------------//
// RNG_Test.hh
// Thomas M. Evans
// Thu May 27 10:45:13 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Some services we will need to test RNG packages.
//---------------------------------------------------------------------------//

#ifndef __rng_test_RNG_Test_hh__
#define __rng_test_RNG_Test_hh__

#include <string>
#include <iostream>

namespace rtt_rng_test
{
 
//===========================================================================//
// FAILURE LIMIT
//===========================================================================//

inline bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    return false;
}

} // end namespace rtt_rng_test

#endif                          // __rng_test_RNG_Test_hh__

//---------------------------------------------------------------------------//
//                              end of rng/test/RNG_Test.hh
//---------------------------------------------------------------------------//
