//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/ds_test.cc
 * \author Thomas M. Evans
 * \date   Wed Nov  7 15:54:59 2001
 * \brief  ds++ testing utilities
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds_test.hh"

namespace rtt_ds_test
{

//===========================================================================//
// PASS/FAILURE
//===========================================================================//

//lint -e755  do not warn about macro functions not being used by test code.
//lint -e759  do not warn about moving function definitions away from header
//            file. 
//lint -e714  do not warn about these functions not being used.
//lint -e552  symbol X not accessed.

bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//

bool fail(int line, const char *file)
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

} // end namespace rtt_ds_test

//---------------------------------------------------------------------------//
//                              end of ds_test.cc
//---------------------------------------------------------------------------//
