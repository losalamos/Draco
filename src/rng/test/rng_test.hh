//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/rng_test.hh
 * \author Thomas M. Evans
 * \date   Mon Dec 17 16:04:59 2001
 * \brief  rng testing services.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __rng_test_hh__
#define __rng_test_hh__

#include <iostream>
#include <string>

namespace rtt_rng_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_rng_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_rng_test::passed will have its default value of true)

// Needless to say, these can be used in many different combinations or
// ways.  We do not constrain draco tests except that the output must be of
// the form "Test: pass/fail"

bool fail(int line);

bool fail(int line, char *file);

bool pass_msg(const std::string &);

bool fail_msg(const std::string &);

//---------------------------------------------------------------------------//
// PASSING CONDITIONALS
//---------------------------------------------------------------------------//

extern bool passed;

} // end namespace rtt_rng_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_rng_test::fail(__LINE__);
#define FAILURE    rtt_rng_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_rng_test::pass_msg(a);
#define FAILMSG(a) rtt_rng_test::fail_msg(a);

#endif                          // __rng_test_hh__

//---------------------------------------------------------------------------//
//                              end of rng/test/rng_test.hh
//---------------------------------------------------------------------------//
