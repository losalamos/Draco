//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstScalarUnitTest.cc
 * \author Kelly Thompson
 * \date   Thu May 18 17:17:24 2006
 * \brief  
 * \note   Copyright 2006 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <map>

#include "../ScalarUnitTest.hh"
#include "../Release.hh"

using namespace std;
using namespace rtt_dsxx;

// Provide old style call to pass/fail macros.  Use object name unitTest for
// this unit test.
#define PASSMSG(a) unitTest.passes(a)
#define ITFAILS    unitTest.failure(__LINE__);
#define FAILURE    unitTest.failure(__LINE__, __FILE__);
#define FAILMSG(a) unitTest.failure(a);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstOne( UnitTest &unitTest )
{
    unitTest.passes("Looks like the passes member function is working.");
    PASSMSG("Looks like the PASSMSG macro is working as a member function.");
    
    return;
}

//---------------------------------------------------------------------------//
// Inclusion of this test will cause failures in the automatic regression
// system because the python filter clues on the keyword "failed".  Only run
// this for stand alone testing
void tstTwo( UnitTest &unitTest )
{
//     cout << "Executing tstTwo:  testing failure functions.\n"
//          << "Output will be saved to a buffer." << endl;
    unitTest.failure("Looks like the failure member function is working.");
    FAILMSG("Looks like the FAILMSG macro is working.");
    ITFAILS;
    FAILURE;

    // Kill report of failures
    unitTest.reset();

    // We need at least one pass.
    PASSMSG("Done with tstTwo.");
    return;
}

//---------------------------------------------------------------------------//
void tstTwoCheck( UnitTest &unitTest, std::ostringstream & msg )
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

//     cout << "The messages from tstTwo contained the following words/occurances."
//           << endl;
//     // print the word_list
//     for( std::map<string,unsigned>::iterator it = word_list.begin();
//          it != word_list.end(); ++it)
//     {
//         cout << it->first << ": " << it->second << endl;
//     }
    
    // Check the list of occurances against the expected values
    if( word_list[ string("Test") ] == 5 )
        unitTest.passes("Found 5 occurances of \"Test\"");
    if( word_list[ string("failed") ] != 4 )
        unitTest.failure("Did not find 4 occurances of failure.");
    if( word_list[ string("FAILMSG") ] != 1 )
        unitTest.failure("Found 1 occurance of \"FAILMSG\"");
    if( word_list[ string("failure") ] != 1 )
        unitTest.failure("Found 1 occurance of \"failure\"");
    if( word_list[ string("macro") ] == 1 )
        unitTest.passes("Found 1 occurance of \"macro\"");
    if( word_list[ string("working") ] == 2 )
        unitTest.passes("Found 2 occurances of \"working\"");
    
    return;
}

//---------------------------------------------------------------------------//
void tstVersion( UnitTest & ut, int & argc, char **& argv )
{
    // build the command that contains "--version"
    string cmd;
    for( int ic=0; ic<argc; ++ic )
        cmd += " " + string( argv[0] );
    cmd += " --version";
    
    system( cmd.c_str() );
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    try
    {
        // Test ctor for ScalarUnitTest (also tests UnitTest ctor and member
        // function setTestName).
        ScalarUnitTest ut( argc, argv, release );
        tstOne(ut);

        // Silent version.
        std::ostringstream messages;
        ScalarUnitTest sut( argc, argv, release, messages );
        tstTwo(sut);

        tstTwoCheck( ut, messages );
        if( argc == 1 )
        {
            // Test --version option.
            tstVersion( ut, argc, argv );
        }

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
//                        end of tstScalarUnitTest.cc
//---------------------------------------------------------------------------//
