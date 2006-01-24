//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   sf/test/sf_test.hh
 * \author Thomas M. Evans
 * \date   Tue Mar 26 12:34:44 2002
 * \brief  sf test infrastructure.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __sf_test_hh__
#define __sf_test_hh__

#include <iostream>
#include <string>

namespace rtt_sf_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_sf_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_sf_test::passed will have its default value of true)

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

} // end namespace rtt_sf_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_sf_test::fail(__LINE__);
#define FAILURE    rtt_sf_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_sf_test::pass_msg(a);
#define FAILMSG(a) rtt_sf_test::fail_msg(a);

#endif                          // __sf_test_hh__

//---------------------------------------------------------------------------//
//                              end of sf/test/sf_test.hh
//---------------------------------------------------------------------------//
