//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diffusion/Tester.cc
 * \author Randy M. Roberts
 * \date   Mon Feb  7 10:00:15 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Tester.hh"

#include <iostream>
#include <sstream>
#include <string>

namespace rtt_diffusion
{

void Tester::testAssert(bool passed, const std::string &msg)
{
    if (!passed)
    {
	os() << Name() << " failed: " << msg << std::endl;
	passed_m = false;
    }
}

void Tester::testAssert(bool passed, const std::string &msg,
			const std::string &file, const std::string &line)
{
    testAssert(passed, msg + " in " + file + " line: " + line);
}

void Tester::testAssert(bool passed, const std::string &msg,
			const std::string &file, int line)
{
    std::ostringstream strline;
    strline << line;
    testAssert(passed, msg, file, strline.str());
}

void Tester::run()
{
    if (checkForVersion())
	return;
    
    try
    {
	os() << "Initiating " << Name() << " test"
		  << std::endl;
	runTest();
    }
    catch( rtt_dsxx::assertion& a )
    {
	os() << "Failed assertion: " << a.what() << std::endl;
	setPassed(false);
    }

    // Print the status of the test.

    printStatus(Name());
}

} // end namespace rtt_diffusion

//---------------------------------------------------------------------------//
//                              end of Tester.cc
//---------------------------------------------------------------------------//
