//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UnitTestFrame/mainDriver.cc
 * \author Randy M. Roberts
 * \date   Fri Feb 25 10:12:53 2000
 * \brief  Main driver for the UnitTestFrame
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
    C4::Init(argc, argv);

    std::string status;

    using rtt_UnitTestFrame::TestApp;

    try
    {
	TestApp &theTestApp = TestApp::theTestApp(argc, argv);
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

    TestApp &theTestApp = TestApp::theTestApp();
    status = theTestApp.run();
    
    std::copy(theTestApp.messageList().begin(),
	      theTestApp.messageList().end(),
	      std::ostream_iterator<std::string>(theTestApp.os(), "\n"));

    if (status.size() > 0)
	theTestApp.os() << status << std::endl;

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of mainDriver.cc
//---------------------------------------------------------------------------//
