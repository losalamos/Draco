//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <pkg>/<spkg>.hh
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef <spkg>_test_hh
#define <spkg>_test_hh

#include <iostream>
#include <string>

namespace <namespace>
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set <namespace>::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (<namespace>::passed will have its default value of true)

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

} // end namespace <namespace>

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS      <namespace>::fail(__LINE__);
#define FAILURE      <namespace>::fail(__LINE__, __FILE__);
#define PASSMSG(a)   <namespace>::pass_msg(a);
#define FAILMSG(a)   <namespace>::fail_msg(a);
#define UNIT_TEST(x) <namespace>::unit_test(x, __LINE__, __FILE__)
    
#endif // <spkg>_test_hh

//---------------------------------------------------------------------------//
//     end of <pkg>/<spkg>.hh
//---------------------------------------------------------------------------//
