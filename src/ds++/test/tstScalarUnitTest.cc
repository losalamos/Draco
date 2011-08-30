//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstScalarUnitTest.cc
 * \author Kelly Thompson
 * \date   Thu May 18 17:17:24 2006
 * \brief  Unit test for the ds++ classes UnitTest and ScalarUnitTest.
 * \note   Copyright © 2006-2010 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include <iostream>
#include <sstream>
#include <map>
#include <cstdlib>
#include <fstream>

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

void tstTwo( UnitTest &unitTest )
{
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
void tstTwoCheck( UnitTest &unitTest, ostringstream & msg )
{
    bool verbose(true);
    map<string,unsigned> word_list(
        UnitTest::get_word_count( msg, verbose ) );

    // Check the list of occurances against the expected values
    if( word_list[ string("Test") ] == 6 )
        unitTest.passes("Found 6 occurances of \"Test\"");
    else
        unitTest.failure("Did not find expected number of occurances of \"Test\"");

    if( word_list[ string("failed") ] != 4 )
        unitTest.failure("Did not find 4 occurances of failure.");
    if( word_list[ string("FAILMSG") ] != 1 )
        unitTest.failure("Found 1 occurance of \"FAILMSG\"");
    if( word_list[ string("failure") ] != 1 )
        unitTest.failure("Found 1 occurance of \"failure\"");

    if( word_list[ string("macro") ] == 1 )
        unitTest.passes("Found 1 occurance of \"macro\"");
    else
        unitTest.failure("Did not find expected number of occurances of \"macro\"");

    if( word_list[ string("working") ] == 2 )
        unitTest.passes("Found 2 occurances of \"working\"");
    else
        unitTest.failure("Did not find expected number of occurances of \"working\"");
    
    return;
}

//---------------------------------------------------------------------------//
void tstVersion( UnitTest & /* ut */, int & argc, char **& argv )
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
void tstGetWordCountFile( UnitTest & unitTest )
{
    cout << "\ntstGetWordCountFile...\n" << endl;
    
    // Generate a text file
    string filename("tstScalarUnitTest.sample.txt");
    ofstream myfile( filename.c_str() );
    if( myfile.is_open() )
    {
        myfile << "This is a text file.\n"
               << "Used by tstScalarUnitTest::tstGetWordCountFile\n\n"
               << "foo bar baz\n"
               << "foo bar\n"
               << "foo\n\n";
        myfile.close();            
    }

    // Now read the file and parse the contents:
    map< string, unsigned > word_list(
        UnitTest::get_word_count( filename, false ) );

    // Some output
    cout << "The world_list has the following statistics (word, count):\n"
         << endl;
    for( map<string,unsigned>::iterator it=word_list.begin();
         it != word_list.end(); ++it )
        cout << it->first << "\t::\t" << it->second << endl;

    // Spot checks on file contents:
    if( word_list[ string("This") ] != 1 )               ITFAILS;
    if( word_list[ string("foo" ) ] != 3 )               ITFAILS;
    if( word_list[ string("bar" ) ] != 2 )               ITFAILS;
    if( word_list[ string("baz" ) ] != 1 )               ITFAILS;

    if( unitTest.numFails == 0 )
    {
        cout << endl;
        PASSMSG("Successfully parsed text file and generated word_list");
    }
    return;
}

//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    try
    {
        // Test ctor for ScalarUnitTest (also tests UnitTest ctor and member
        // function setTestName).
        ScalarUnitTest ut( argc, argv, release );
        tstOne(ut);

        // Silent version.
        ostringstream messages;
        ScalarUnitTest sut( argc, argv, release, messages );
        tstTwo(sut);

        tstTwoCheck( ut, messages );

        tstGetWordCountFile( ut );
        
        if( argc == 1 )
        {
            // Test --version option.
            tstVersion( ut, argc, argv );
        }

    }
    catch( rtt_dsxx::assertion &err )
    {
        string msg = err.what();
        if( msg != string( "Success" ) )
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
