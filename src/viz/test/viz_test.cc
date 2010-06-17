//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/test/viz_test.cc
 * \author Rob Lowrie
 * \date   Mon Mar  8 10:29:48 2004
 * \brief  Test stuff for viz
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include "viz_test.hh"

namespace rtt_viz_test
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

bool fail(int line, char const * file)
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

void unit_test(const bool pass, int line, char const * file)
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

} // end namespace rtt_viz_test

//---------------------------------------------------------------------------//
//               end of viz_test.cc
//---------------------------------------------------------------------------//
