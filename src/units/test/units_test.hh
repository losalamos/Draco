//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   units/test/units_test.hh
 * \author U-LUMEN\kellyt
 * \date   Wed Oct  8 13:49:47 2003
 * \brief  
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_units_test_test_hh
#define rtt_units_test_test_hh

#include "ds++/config.h"
#include <iostream>
#include <string>

namespace rtt_units_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_units_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_units_test::passed will have its default value of true)

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

DLL_PUBLIC extern bool passed;

} // end namespace rtt_units_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_units_test::fail(__LINE__);
#define FAILURE    rtt_units_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_units_test::pass_msg(a);
#define FAILMSG(a) rtt_units_test::fail_msg(a);

#endif // rtt_units_test_test_hh

//---------------------------------------------------------------------------//
//     end of units/test/units_test.hh
//---------------------------------------------------------------------------//
