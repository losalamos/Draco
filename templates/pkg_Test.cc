//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <pkg>/<spkg>.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include "<spkg>.hh"

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

void unit_test(const bool pass, int line, char *file)
{
    if ( pass )
    {
        std::cout << "Test: passed\n";
    }
    else
    {
	fail(line, file);
    }
}

//---------------------------------------------------------------------------//
// BOOLEAN PASS FLAG
//---------------------------------------------------------------------------//

bool passed = true;

} // end namespace rtt_<spkg>

//---------------------------------------------------------------------------//
//               end of <spkg>.cc
//---------------------------------------------------------------------------//
