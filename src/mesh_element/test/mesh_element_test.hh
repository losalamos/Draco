//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh_element/test/mesh_element_test.hh
 * \author Kelly Thompson
 * \date   Mon May 24 16:58:08 2004
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef mesh_element_test_test_hh
#define mesh_element_test_test_hh

#include <iostream>
#include <string>

namespace rtt_mesh_element_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_mesh_element_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_mesh_element_test::passed will have its default value of true)

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

} // end namespace rtt_mesh_element_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS      rtt_mesh_element_test::fail(__LINE__);
#define FAILURE      rtt_mesh_element_test::fail(__LINE__, __FILE__);
#define PASSMSG(a)   rtt_mesh_element_test::pass_msg(a);
#define FAILMSG(a)   rtt_mesh_element_test::fail_msg(a);
#define UNIT_TEST(x) rtt_mesh_element_test::unit_test(x, __LINE__, __FILE__)
    
#endif // mesh_element_test_test_hh

//---------------------------------------------------------------------------//
//     end of mesh_element/test/mesh_element_test.hh
//---------------------------------------------------------------------------//
