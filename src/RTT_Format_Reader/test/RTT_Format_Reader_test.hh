//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/test/RTT_Format_Reader_test.hh
 * \author Thomas M. Evans
 * \date   Tue Mar 26 17:12:55 2002
 * \brief  RTT_Format_Reader testing infrastructure.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_Reader_test_hh__
#define __RTT_Format_Reader_test_hh__

#include <iostream>
#include <string>

namespace rtt_RTT_Format_Reader_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_RTT_Format_Reader_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_RTT_Format_Reader_test::passed will have its default value of true)

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

} // end namespace rtt_RTT_Format_Reader_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_RTT_Format_Reader_test::fail(__LINE__);
#define FAILURE    rtt_RTT_Format_Reader_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_RTT_Format_Reader_test::pass_msg(a);
#define FAILMSG(a) rtt_RTT_Format_Reader_test::fail_msg(a);

#endif                          // __RTT_Format_Reader_test_hh__

//---------------------------------------------------------------------------//
//                              end of RTT_Format_Reader/test/RTT_Format_Reader_test.hh
//---------------------------------------------------------------------------//
