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

#include "TestAppBase.hh"

int main(int argc, char *argv[])
{
    using rtt_UnitTestFrame::TestAppBase;

    // This does the whole enchilada!
    // It should initialize the process, create the TestApp,
    // run the TestApp, and finalize the process.
    
    return TestAppBase::driveTest(argc, argv);
}

//---------------------------------------------------------------------------//
//                              end of mainDriver.cc
//---------------------------------------------------------------------------//
