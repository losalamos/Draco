//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshReaders_Services/test/meshReaders_Services_test.cc
 * \author Thomas M. Evans
 * \date   Wed Mar 27 14:53:05 2002
 * \brief  meshReaders_Services testing infrastructure.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "meshReaders_Services_test.hh"
#include <iostream>

namespace rtt_meshReaders_Services_test
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

} // end namespace rtt_meshReaders_Services_test

//---------------------------------------------------------------------------//
//                              end of meshReaders_Services_test.cc
//---------------------------------------------------------------------------//
