//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   norms/test/norms_test.hh
 * \author Rob Lowrie
 * \date   Fri Jan 14 09:10:18 2005
 * \brief  header for utilities.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_norms_test_hh
#define rtt_norms_test_hh

#include <iostream>
#include <string>

namespace rtt_norms_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_norms_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_norms_test::passed will have its default value of true)

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

} // end namespace rtt_norms_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS      rtt_norms_test::fail(__LINE__);
#define FAILURE      rtt_norms_test::fail(__LINE__, __FILE__);
#define PASSMSG(a)   rtt_norms_test::pass_msg(a);
#define FAILMSG(a)   rtt_norms_test::fail_msg(a);
#define UNIT_TEST(x) rtt_norms_test::unit_test(x, __LINE__, __FILE__)
    
#endif // rtt_norms_test_hh

//---------------------------------------------------------------------------//
//     end of norms/test/norms_test.hh
//---------------------------------------------------------------------------//
