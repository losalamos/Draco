//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/mainDriver.cc
 * \author Randy M. Roberts
 * \date   Fri Feb 25 10:12:53 2000
 * \brief  Main driver for the UnitTestFrame
 *
 * This file is used to provide the main program that drives unit tests
 * within the unit test framework.  The user of this framework need only
 * derive a class from TestApp, supply any pure virtual methods specified
 * by TestApp, provide a definition for TestApp::create(), and link in
 * this library.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "TestApp.hh"

#include "c4/global.hh"

#include <exception>
#include <iostream>

int main(int argc, char *argv[])
{
    // This is the call to kick off C4.
    
    C4::Init(argc, argv);

    // Decide wether we want to use try blocks or cause core
    // dumps. The --core flag on the command line determines that core dumps,
    // instead of try blocks will be used.
    
    bool coreDump = false;
    
    for (int i=1; i<argc; ++i)
	if (std::string(argv[i]) == "--core")
	    coreDump = true;

    std::string status;

    using rtt_UnitTestFrame::TestApp;

    if (coreDump)
    {
	try
	{
	    // Create and access the TestApp derived class singleton.
	    
	    TestApp &theTestApp = TestApp::theTestApp(argc, argv, std::cout);
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
	// Create and access the TestApp derived class singleton.
	
	TestApp &theTestApp = TestApp::theTestApp(argc, argv, std::cout);
    }

    // This access of the TestApp derived class object does **not** create
    // a new one... it's a singleton for goodness sake!
    
    TestApp &theTestApp = TestApp::theTestApp();

    // Let's run the unit test.
    // The base class run method will check for "--core" and decide whether
    // to run the test within a try block.
    
    status = theTestApp.run();

    // Copy all of the generated pass/fail messages to the pass/fail ostream.
    
    std::copy(theTestApp.messageList().begin(),
	      theTestApp.messageList().end(),
	      std::ostream_iterator<std::string>(theTestApp.os(), "\n"));

    // If the run method returned a string let's copy it to the pass/fail
    // ostream.
    
    if (status.size() > 0)
	theTestApp.os() << status << std::endl;

    // This is the call to kick out C4.
    
    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of mainDriver.cc
//---------------------------------------------------------------------------//
