//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders/test/meshReaders_test.hh
 * \author Thomas M. Evans
 * \date   Tue Mar 26 16:05:39 2002
 * \brief  meshReaders testing infrastructure.
 * \note   Copyright (C) 2002-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_test_hh__
#define __meshReaders_test_hh__

#include "ds++/config.h"
#include <iostream>
#include <string>

namespace rtt_meshReaders_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_meshReaders_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_meshReaders_test::passed will have its default value of true)

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

} // end namespace rtt_meshReaders_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_meshReaders_test::fail(__LINE__);
#define FAILURE    rtt_meshReaders_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_meshReaders_test::pass_msg(a);
#define FAILMSG(a) rtt_meshReaders_test::fail_msg(a);

#endif // __meshReaders_test_hh__

//---------------------------------------------------------------------------//
// end of meshReaders/test/meshReaders_test.hh
//---------------------------------------------------------------------------//
