//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/TestApp.cc
 * \author Randy M. Roberts
 * \date   Fri Feb 25 07:34:16 2000
 * \brief  The abstract base from which unit test applications can be derived.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestApp.hh"
#include "PassFailStream.hh"
#include <algorithm>

using std::string;

namespace rtt_UnitTestFrame
{

rtt_dsxx::SP<TestApp> TestApp::theSPTestApp;

TestApp::TestApp(int &argc,  char *argv[], std::ostream &os_in)
    : os_m(os_in), passed_m(true)
{
    // Create the argument list.
    
    for (int i=1; i<argc; ++i)
	argList_m.push_back(argv[i]);
}

string TestApp::run()
{
    // If the user specified --version on the command line,
    // we want to return the version... That's it!
    
    if (std::find(argList_m.begin(), argList_m.end(), "--version")
	!= argList_m.end())
    {
	return name() + " version: " + version();
    }

    // If the user specified --core on the command line,
    // we do not want to put runTest() call within a try block.
    
    if (std::find(argList_m.begin(), argList_m.end(), "--core")
	!= argList_m.end())
    {
	return runTest();
    }

    // If the user did not specify --core then we must
    // call runTest() within a try block.
    
    try
    {
	return runTest();
    }
    catch(rtt_dsxx::assertion& a)
    {
	fail() << a.what();
    }
    catch(std::exception& a)
    {
	fail() <<"std::exception: " << a.what();
    }
    return string("Test ") + name() + " incomplete.";
}

PassFailStream TestApp::pass(const string &str)
{
    // Return a passing PassFailStream.
    
    return PassFailStream(*this, str, true);			
}

PassFailStream TestApp::fail(const string &str)
{
    // Return a failing PassFailStream.
    
    return PassFailStream(*this, str, false);			
}

string TestApp::getNextArg(const string &arg) const
{
    // Get the argument after the argument determined by arg.
    
    std::list<string>::const_iterator argit;
    argit = std::find(argList_m.begin(), argList_m.end(), arg);

    // If we find the argument, and its not the last argument...
    
    if (argit != argList_m.end() && ++argit != argList_m.end())
	return *argit;

    return "";
}

} // end namespace rtt_UnitTestFrame


//---------------------------------------------------------------------------//
//                              end of TestApp.cc
//---------------------------------------------------------------------------//
