//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/TestAppBase.cc
 * \author Randy M. Roberts
 * \date   Fri Feb 25 07:34:16 2000
 * \brief  The abstract base from which unit test applications can be derived.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestAppBase.hh"
#include "PassFailStream.hh"
#include "ProcessPolicy.hh"

#include <algorithm>
#include <exception>
#include <iostream>

using std::string;

namespace rtt_UnitTestFrame
{

rtt_dsxx::SP<TestAppBase> TestAppBase::theSPTestAppBase;

TestAppBase::TestAppBase(int &argc,  char *argv[], std::ostream &os_in)
    : os_m(os_in), passed_m(true)
{
    // Create the argument list.
    
    for (int i=1; i<argc; ++i)
	argList_m.push_back(argv[i]);
}

int TestAppBase::driveTest(int &argc, char *argv[])
{
    // This is the call to kick off C4.
    
    ProcessInit(argc, argv);

    // Decide wether we want to use try blocks or cause core
    // dumps. The --core flag on the command line determines that core dumps,
    // instead of try blocks will be used.
    
    bool coreDump = false;
    
    for (int i=1; i<argc; ++i)
	if (std::string(argv[i]) == "--core")
	    coreDump = true;

    std::string status;

    using rtt_UnitTestFrame::TestAppBase;

    if (coreDump)
    {
	try
	{
	    // Create and access the TestAppBase derived class singleton.
	    
	    TestAppBase &theTestAppBase =
                TestAppBase::theTestAppBase(argc, argv, std::cout);
	}
	catch(rtt_dsxx::assertion& a)
	{
	    status = std::string("Caught Assertion in Test Construction: ")
		+ a.what();
	}
	catch(std::exception& a)
	{
	    status = std::string("Caught Exception in Test Construction: ")
		+ a.what();
	}
    }
    else
    {
	// Create and access the TestAppBase derived class singleton.
	
	TestAppBase &theTestAppBase = TestAppBase::theTestAppBase(argc, argv,
                                                                  std::cout);
    }

    // This access of the TestAppBase derived class object does **not** create
    // a new one... it's a singleton for goodness sake!
    
    TestAppBase &theTestAppBase = TestAppBase::theTestAppBase();

    // Let's run the unit test.
    // The base class run method will check for "--core" and decide whether
    // to run the test within a try block.
    
    status = theTestAppBase.run();

    // Copy all of the generated pass/fail messages to the pass/fail ostream.
    
    std::copy(theTestAppBase.messageList().begin(),
	      theTestAppBase.messageList().end(),
	      std::ostream_iterator<std::string>(theTestAppBase.os(), "\n"));

    // If the run method returned a string let's copy it to the pass/fail
    // ostream.
    
    if (status.size() > 0)
	theTestAppBase.os() << status << std::endl;

    ProcessFinalize();
    
    // return theTestAppBase.passed() ? 0 : 1;
    return 0;
}

string TestAppBase::run()
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

PassFailStream TestAppBase::pass(const string &str)
{
    // Return a passing PassFailStream.
    
    return PassFailStream(*this, str, true);			
}

PassFailStream TestAppBase::fail(const string &str)
{
    // Return a failing PassFailStream.
    
    return PassFailStream(*this, str, false);			
}

string TestAppBase::getNextArg(const string &arg) const
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
//                              end of TestAppBase.cc
//---------------------------------------------------------------------------//
