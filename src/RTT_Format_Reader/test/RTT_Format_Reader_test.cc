//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   RTT_Format_Reader/test/RTT_Format_Reader_test.cc
 * \author Thomas M. Evans
 * \date   Tue Mar 26 17:12:55 2002
 * \brief  RTT_Format_Reader testing infrastructure.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "RTT_Format_Reader_test.hh"
#include <iostream>

namespace rtt_RTT_Format_Reader_test
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

} // end namespace rtt_RTT_Format_Reader_test

//---------------------------------------------------------------------------//
//                              end of RTT_Format_Reader_test.cc
//---------------------------------------------------------------------------//
