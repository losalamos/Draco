//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/Tester.cc
 * \author Randy M. Roberts
 * \date   Wed Nov 24 13:09:05 1999
 * \brief  Implementation for Tester class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Tester.hh"

#include <iostream>
#include <sstream>
#include <string>

namespace rtt_meshTest
{

void Tester::testassert(bool passed, const std::string &msg)
{
    if (!passed)
    {
	os_m << name_m << " failed: " << msg << std::endl;
	passed_m = false;
    }
}

void Tester::testassert(bool passed, const std::string &msg,
			const std::string &file, const std::string &line)
{
    testassert(passed, msg + " in " + file + " line: " + line);
}

void Tester::testassert(bool passed, const std::string &msg,
			const std::string &file, int line)
{
    std::ostringstream strline;
    strline << line;
    testassert(passed, msg, file, strline.str());
}

} // end namespace rtt_meshTest

//---------------------------------------------------------------------------//
//                              end of Tester.cc
//---------------------------------------------------------------------------//
