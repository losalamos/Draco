//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/c4_test.hh
 * \author Thomas M. Evans
 * \date   Mon Mar 25 15:30:56 2002
 * \brief  c4 package testing functions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __c4_test_hh__
#define __c4_test_hh__

#include "ds++/config.h"
#include <iostream>
#include <string>

namespace rtt_c4_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_c4_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_c4_test::passed will have its default value of true)

// Needless to say, these can be used in many different combinations or
// ways.  We do not constrain draco tests except that the output must be of
// the form "Test: pass/fail"

DLL_PUBLIC bool fail(int line);

DLL_PUBLIC bool fail(int line, char *file);

DLL_PUBLIC bool pass_msg(const std::string &);

DLL_PUBLIC bool fail_msg(const std::string &);

//---------------------------------------------------------------------------//
// PASSING CONDITIONALS
//---------------------------------------------------------------------------//

extern DLL_PUBLIC bool passed;

} // end namespace rtt_c4_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_c4_test::fail(__LINE__);
#define FAILURE    rtt_c4_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_c4_test::pass_msg(a);
#define FAILMSG(a) rtt_c4_test::fail_msg(a);

#endif                          // __c4_test_hh__

//---------------------------------------------------------------------------//
//                              end of c4/test/c4_test.hh
//---------------------------------------------------------------------------//
