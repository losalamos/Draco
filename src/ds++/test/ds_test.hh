//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/ds_test.hh
 * \author Thomas M. Evans
 * \date   Wed Nov  7 15:54:59 2001
 * \brief  ds++ testing utilities.
 * \note   Copyright (c) 1997-2010 Los Alamos National Security, LLC 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __ds_test_hh__
#define __ds_test_hh__

#include "ds++/config.h"
#include <iostream>
#include <string>

namespace rtt_ds_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_ds_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_ds_test::passed will have its default value of true)

// Needless to say, these can be used in many different combinations or
// ways.  We do not constrain draco tests except that the output must be of
// the form "Test: pass/fail"

DLL_PUBLIC bool fail(int line);

DLL_PUBLIC bool fail(int line, char const *file);

DLL_PUBLIC bool pass_msg(const std::string &);

DLL_PUBLIC bool fail_msg(const std::string &);

DLL_PUBLIC void unit_test(const bool pass, int line, char const *file);

//---------------------------------------------------------------------------//
// PASSING CONDITIONALS
//---------------------------------------------------------------------------//

extern 
DLL_PUBLIC bool passed;

} // end namespace rtt_ds_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_ds_test::fail(__LINE__)
#define FAILURE    rtt_ds_test::fail(__LINE__, __FILE__)
#define PASSMSG(a) rtt_ds_test::pass_msg(a)
#define FAILMSG(a) rtt_ds_test::fail_msg(a)
#define UNIT_TEST(x) rtt_ds_test::unit_test(x, __LINE__, __FILE__)

#endif                          // __ds_test_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/test/ds_test.hh
//---------------------------------------------------------------------------//
