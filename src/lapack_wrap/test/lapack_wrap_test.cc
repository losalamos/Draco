//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   lapack_wrap/test/lapack_wrap_test.cc
 * \author Thomas M. Evans
 * \date   Thu Aug 29 11:30:25 2002
 * \brief  lapack_wrap testing services.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "lapack_wrap_test.hh"
#include <iostream>

namespace rtt_lapack_wrap_test
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

} // end namespace rtt_lapack_wrap_test

//---------------------------------------------------------------------------//
//                              end of lapack_wrap_test.cc
//---------------------------------------------------------------------------//
