//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/test/CDI_Test.cc
 * \author Thomas M. Evans
 * \date   Thu Sep  6 18:12:20 2001
 * \brief  Testing services for cdi_analytic
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "CDI_Test.hh"
#include <iostream>

namespace rtt_cdi_analytic_test
{

//===========================================================================//
// FAILURE LIMITS
//===========================================================================//

bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    return false;
}

bool fail(int line, char *file)
{
    std::cout << "Test: failed on line " << line << " in " << file
	      << std::endl;
    return false;
}

} // end namespace rtt_cdi_analytic_test

//---------------------------------------------------------------------------//
//                              end of CDI_Test.cc
//---------------------------------------------------------------------------//
