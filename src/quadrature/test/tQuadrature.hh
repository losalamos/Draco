//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Quadrature/tQuadrature.hh
 * \author Kelly Thompson
 * \date   Mon Mar 6 13:41:03 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __Quadrature_test_tQuadrature_hh__
#define __Quadrature_test_tQuadrature_hh__

#include <string>

namespace rtt_quadrature_test
{
 
using std::string;

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

} // end namespace rtt_Quadrature_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_quadrature_test::fail(__LINE__);
#define FAILURE    rtt_quadrature_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_quadrature_test::pass_msg(a);
#define FAILMSG(a) rtt_quadrature_test::fail_msg(a);

#endif // __Quadrature_test_tQuadrature_hh__

//---------------------------------------------------------------------------//
//		      end of Quadrature/tQuadrature.hh
//---------------------------------------------------------------------------//
