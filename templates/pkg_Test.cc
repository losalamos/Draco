//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <pkg>/<spkg>.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "<spkg>.hh"
#include <iostream>

namespace rtt_<spkg>
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

} // end namespace rtt_<spkg>

//---------------------------------------------------------------------------//
//                              end of <spkg>.cc
//---------------------------------------------------------------------------//
