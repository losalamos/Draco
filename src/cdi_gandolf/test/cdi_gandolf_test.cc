//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/test/cdi_gandolf_test.cc
 * \author Thomas M. Evans
 * \date   Fri Oct 12 15:36:36 2001
 * \brief  cdi_gandolf test functions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_gandolf_test.hh"
#include <iostream>

namespace rtt_cdi_gandolf_test
{

//===========================================================================//
// PASS/FAILURE
//===========================================================================//

bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//

bool fail(int line, char *file)
{
    std::cout << "Test: failed on line " << line << " in " << file
	      << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//

bool pass_msg(const std::string &passmsg)
{
    std::cout << "Test: passed" << std::endl;
    std::cout << "     " << passmsg << std::endl;
    return true;
}

//---------------------------------------------------------------------------//

bool fail_msg(const std::string &failmsg)
{
    std::cout << "Test: failed" << std::endl;
    std::cout << "     " << failmsg << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//
// BOOLEAN PASS FLAG
//---------------------------------------------------------------------------//

bool passed = true;

} // end namespace rtt_cdi_gandolf_test

//---------------------------------------------------------------------------//
//                              end of cdi_gandolf_test.cc
//---------------------------------------------------------------------------//
