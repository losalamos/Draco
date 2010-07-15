//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstApplicationUnitTest.cc
 * \author Kelly Thompson
 * \date   Tue Jun  6 15:03:08 2006
 * \brief  Test the Draco class ApplicationUnitTest
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <map>

#include "../Release.hh"
#include "../ApplicationUnitTest.hh"

using namespace std;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstOne( ApplicationUnitTest &unitTest )
{
    string const extraArg( "hello" );
    unitTest.addCommandLineArgument( extraArg );
    cout << ">>> Executing unitTest.runTests()..." << endl;
    unitTest.runTests();

    //! \bug Consider using Boost or other 3rd party library to aid with
    // file path manipulation, including finding the cwd.
    
    string const logFilename( unitTest.logFileName() );
    ostringstream msg;
    msg << "phw_hello-"<< unitTest.nodes() <<".out";
    string const expLogFilename( msg.str() );

    // Find the expected filename (no path) in the real filename
    size_t pos( logFilename.find( expLogFilename ) );
    if( pos != string::npos )
        unitTest.passes( "Found expected log filename (pos != npos)." );
    else
        unitTest.failure( "Did not find expected log filename (pos == npos)." );
    cout << endl;
    return;
}

//---------------------------------------------------------------------------//

void tstTwo( ApplicationUnitTest &unitTest )
{
    // This test is designed to fail.
    
    std::string const extraArg;
    std::cout << ">>> Executing unitTest.runTest( extraArg )..." << std::endl;
    if( unitTest.runTest( extraArg ) )
        unitTest.failure("Found problems when running phw.");
    else
        unitTest.passes("Successfully ran phw.");

    // Kill fail flag (we expected this failure).
    unitTest.reset();
    // We need at least one pass.
    unitTest.passes("Done with tstTwo.");
    if( unitTest.allTestsPass() )
        unitTest.passes("All tests pass.");
    else
        unitTest.failure("Some tests failed.");
    return;
}

//---------------------------------------------------------------------------//

void tstTwoCheck( ApplicationUnitTest &unitTest, std::ostringstream & msg )
{
    using rtt_dsxx::UnitTest;
    
    std::map<string,unsigned> word_list( UnitTest::get_word_count( msg ) );
    
    // Check the list of occurances against the expected values
    if( word_list[ string("Test") ] == 5 )
        unitTest.passes("Found 5 occurances of \"Test\"");
    if( word_list[ string("failed") ] != 3 )
        unitTest.failure("Did not find 3 occurances of failure.");
    if( word_list[ string("passed") ] == 2 )
        unitTest.passes("Found 2 occurances of \"working\"");
    
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        // Test ctor for ApplicationUnitTest 
        ApplicationUnitTest ut( argc, argv, release, "./phw" );
        tstOne(ut);

        // Silent version.
        std::ostringstream messages;
        ApplicationUnitTest sut(
            argc, argv, release, "./phw", std::list<std::string>(), messages );
        tstTwo(sut);
        tstTwoCheck( ut, messages );
        
        ut.status();
    }
    catch( rtt_dsxx::assertion &err )
    {
        std::string msg = err.what();
        if( msg != std::string( "Success" ) )
        { cout << "ERROR: While testing " << argv[0] << ", "
               << err.what() << endl;
            return 1;
        }
        return 0;
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing " << argv[0] << ", "
             << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing " << argv[0] << ", " 
             << "An unknown exception was thrown" << endl;
        return 1;
    }
    return 0;
}   

//---------------------------------------------------------------------------//
//                        end of tstApplicationUnitTest.cc.cc
//---------------------------------------------------------------------------//
