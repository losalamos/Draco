//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/test/diagnostics_test.hh
 * \author Thomas M. Evans
 * \date   Fri Dec  9 11:08:34 2005
 * \brief  diagnostics testing harness.
 * \note   Copyright (C) 2004-2010 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef diagnostics_test_hh
#define diagnostics_test_hh

#include <iostream>
#include <string>

namespace rtt_diagnostics_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_diagnostics_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_diagnostics_test::passed will have its default value of true)

// Needless to say, these can be used in many different combinations or
// ways.  We do not constrain draco tests except that the output must be of
// the form "Test: pass/fail"

bool fail(int line);

bool fail(int line, char *file);

bool pass_msg(const std::string &);

bool fail_msg(const std::string &);

void unit_test(const bool pass, int line, char *file);

//---------------------------------------------------------------------------//
// PASSING CONDITIONALS
//---------------------------------------------------------------------------//

extern bool passed;

} // end namespace rtt_diagnostics_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS      rtt_diagnostics_test::fail(__LINE__);
#define FAILURE      rtt_diagnostics_test::fail(__LINE__, __FILE__);
#define PASSMSG(a)   rtt_diagnostics_test::pass_msg(a);
#define FAILMSG(a)   rtt_diagnostics_test::fail_msg(a);
#define UNIT_TEST(x) rtt_diagnostics_test::unit_test(x, __LINE__, __FILE__)
    
#endif // diagnostics_test_hh

//---------------------------------------------------------------------------//
//     end of diagnostics/diagnostics_test.hh
//---------------------------------------------------------------------------//
