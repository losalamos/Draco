//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/test/cdi_eospac_test.hh
 * \author Kelly Thompson
 * \date   Mon Apr 2 14:15:57 2001
 * \brief  Header file for cdi_eospac_test.cc and tEospac.cc
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_eospac_test_hh__
#define __cdi_eospac_test_hh__

#include <vector>
#include <string>

namespace rtt_cdi_eospac_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_cdi_eospac_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_cdi_eospac_test::passed will have its default value of true)

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
// DATA EQUIVALENCE FUNCTIONS USED FOR TESTING
//---------------------------------------------------------------------------//

bool match( const double computedValue, const double referenceValue );

//---------------------------------------------------------------------------//

bool match(const std::vector< double >& computedValue, 
	   const std::vector< double >& referenceValue );

//---------------------------------------------------------------------------//

bool match(const std::vector< std::vector<double> >& computedValue, 
	   const std::vector< std::vector<double> >& referenceValue );

} // end namespace rtt_cdi_eospac_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_cdi_eospac_test::fail(__LINE__);
#define FAILURE    rtt_cdi_eospac_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_cdi_eospac_test::pass_msg(a);
#define FAILMSG(a) rtt_cdi_eospac_test::fail_msg(a);

#endif // __cdi_eospac_test_hh__

//---------------------------------------------------------------------------//
// end of cdi_eospac/test/tEospac.hh
//---------------------------------------------------------------------------//

