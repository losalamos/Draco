//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   norms/test/norms_test.cc
 * \author Rob Lowrie
 * \date   Fri Jan 14 09:10:19 2005
 * \brief  utilities for tests.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include "norms_test.hh"

namespace rtt_norms_test
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
    if ( ! pass ) fail(line, file);
}

//---------------------------------------------------------------------------//
// BOOLEAN PASS FLAG
//---------------------------------------------------------------------------//

bool passed = true;

} // end namespace rtt_norms_test

//---------------------------------------------------------------------------//
//               end of norms_test.cc
//---------------------------------------------------------------------------//
