//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/cdi_test.hh
 * \author Thomas M. Evans
 * \date   Tue Oct  9 10:51:39 2001
 * \brief  CDI Test help function prototypes.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_test_hh__
#define __cdi_test_hh__

#include <iostream>
#include <vector>

namespace rtt_cdi_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_cdi_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_cdi_test::passed will have its default value of true)

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

//---------------------------------------------------------------------------//
// CHECK COMPUTED VERSUS EXPECTED VALUES
//---------------------------------------------------------------------------//

bool match(double computedValue, double referenceValue);


bool match(const std::vector< double > &computedValue, 
	   const std::vector< double > &referenceValue );

bool match(const std::vector< std::vector< double > >& computedValue, 
	   const std::vector< std::vector< double > >& referenceValue ); 

} // end namespace rtt_cdi_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_cdi_test::fail(__LINE__);
#define FAILURE    rtt_cdi_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_cdi_test::pass_msg(a);
#define FAILMSG(a) rtt_cdi_test::fail_msg(a);

#endif                          // __cdi_test_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/test/cdi_test.hh
//---------------------------------------------------------------------------//
