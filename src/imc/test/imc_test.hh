//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/test/imc_test.hh
 * \author Thomas M. Evans
 * \date   Wed Nov 14 16:50:26 2001
 * \brief  IMC package testing functions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __imc_test_hh__
#define __imc_test_hh__

#include <iostream>
#include <string>

namespace rtt_imc_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_imc_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_imc_test::passed will have its default value of true)

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

} // end namespace rtt_imc_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_imc_test::fail(__LINE__);
#define FAILURE    rtt_imc_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_imc_test::pass_msg(a);
#define FAILMSG(a) rtt_imc_test::fail_msg(a);

#endif                          // __imc_test_hh__

//---------------------------------------------------------------------------//
//                              end of imc/test/imc_test.hh
//---------------------------------------------------------------------------//
