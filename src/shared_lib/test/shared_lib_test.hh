//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   shared_lib/test/shared_lib_test.hh
 * \author Thomas M. Evans
 * \date   Wed Apr 21 14:31:07 2004
 * \brief  shared_lib testing infrastructure.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved.
 */
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

#ifndef shared_lib_test_test_hh
#define shared_lib_test_test_hh

#include <iostream>
#include <string>

namespace rtt_shared_lib_test {

//============================================================================//
// PASS/FAILURE LIMIT
//============================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_shared_lib_test::passed to false

// These can be used in any combination in a test to print output messages
// if no fail functions are called then the test will pass
// (rtt_shared_lib_test::passed will have its default value of true)

// Needless to say, these can be used in many different combinations or
// ways.  We do not constrain draco tests except that the output must be of
// the form "Test: pass/fail"

bool fail(int line);

bool fail(int line, char *file);

bool pass_msg(const std::string &);

bool fail_msg(const std::string &);

void unit_test(const bool pass, int line, char *file);

//----------------------------------------------------------------------------//
// PASSING CONDITIONALS
//----------------------------------------------------------------------------//

extern bool passed;

} // end namespace rtt_shared_lib_test

//============================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//============================================================================//

#define ITFAILS rtt_shared_lib_test::fail(__LINE__);
#define FAILURE rtt_shared_lib_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_shared_lib_test::pass_msg(a);
#define FAILMSG(a) rtt_shared_lib_test::fail_msg(a);
#define UNIT_TEST(x) rtt_shared_lib_test::unit_test(x, __LINE__, __FILE__)

#endif // shared_lib_test_test_hh

//----------------------------------------------------------------------------//
// end of shared_lib/test/shared_lib_test.hh
//----------------------------------------------------------------------------//
