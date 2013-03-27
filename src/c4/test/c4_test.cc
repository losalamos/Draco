//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/c4_test.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 15:39:13 2002
 * \brief  c4 package test infrastructure.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "c4_test.hh"
#include <iostream>

namespace rtt_c4_test
{

//===========================================================================//
// PASS/FAILURE
//===========================================================================//

DLL_PUBLIC bool fail(int line)
{
    std::cout << "Test: failed on line " << line << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//

DLL_PUBLIC bool fail(int line, char *file)
{
    std::cout << "Test: failed on line " << line << " in " << file
	      << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//

DLL_PUBLIC bool pass_msg(const std::string &passmsg)
{
    std::cout << "Test: passed" << std::endl;
    std::cout << "     " << passmsg << std::endl;
    return true;
}

//---------------------------------------------------------------------------//

DLL_PUBLIC bool fail_msg(const std::string &failmsg)
{
    std::cout << "Test: failed" << std::endl;
    std::cout << "     " << failmsg << std::endl;
    passed = false;
    return false;
}

//---------------------------------------------------------------------------//
// BOOLEAN PASS FLAG
//---------------------------------------------------------------------------//

DLL_PUBLIC bool passed = true;

} // end namespace rtt_c4_test

//---------------------------------------------------------------------------//
//                              end of c4_test.cc
//---------------------------------------------------------------------------//
