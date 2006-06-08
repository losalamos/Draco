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
    std::string const extraArg( "hello" );
    unitTest.addCommandLineArgument( extraArg );
    std::cout << ">>> Executing unitTest.runTests()..." << std::endl;
    unitTest.runTests();

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
    std::map<string,unsigned> word_list;
    string msgbuf( msg.str() );
    string delims(" \n\t:,.;");

    { // Build a list of words found in msgbuf.  Count the number of
      // occurances.
        
        // Find the beginning of the first word.
        string::size_type begIdx = msgbuf.find_first_not_of(delims);
        string::size_type endIdx;
        
        // While beginning of a word found
        while( begIdx != string::npos )
        {
            // search end of actual word
            endIdx = msgbuf.find_first_of( delims, begIdx );
            if( endIdx == string::npos)
                endIdx = msgbuf.length();
            
            // the word is we found is...
            string word( msgbuf, begIdx, endIdx-begIdx );
            
            // add it to the map
            word_list[ word ]++;
            
            // search to the beginning of the next word
            begIdx = msgbuf.find_first_not_of( delims, endIdx );        
        }
    }

//  {
//     cout << "The messages from tstTwo contained the following words/occurances."
//           << endl;
//     // print the word_list
//     for( std::map<string,unsigned>::iterator it = word_list.begin();
//          it != word_list.end(); ++it)
//     {
//         cout << it->first << ": " << it->second << endl;
//     }
//  }
    
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
