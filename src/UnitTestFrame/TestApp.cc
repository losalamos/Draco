//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/TestApp.cc
 * \author Randy M. Roberts
 * \date   Fri Feb 25 07:34:16 2000
 * \brief  
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
    for (int i=1; i<argc; ++i)
	argList_m.push_back(argv[i]);
}

string TestApp::run()
{
    if (std::find(argList_m.begin(), argList_m.end(), "--version")
	!= argList_m.end())
    {
	return name() + " version: " + version();
    }

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
    return PassFailStream(*this, str, true);			
}

PassFailStream TestApp::fail(const string &str)
{
    return PassFailStream(*this, str, false);			
}

string TestApp::getNextArg(const string &arg) const
{
    std::list<string>::const_iterator argit;
    argit = std::find(argList_m.begin(), argList_m.end(), arg);

    if (argit != argList_m.end() && ++argit != argList_m.end())
	return *argit;

    return "";
}

} // end namespace rtt_UnitTestFrame


//---------------------------------------------------------------------------//
//                              end of TestApp.cc
//---------------------------------------------------------------------------//
