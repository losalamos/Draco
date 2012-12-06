//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstScalarUnitTest.cc
 * \author Kelly Thompson
 * \date   Thu May 18 17:17:24 2006
 * \brief  Unit test for the ds++ classes UnitTest and ScalarUnitTest.
 * \note   Copyright (C) 2006-2012 Los Alamos National Security, LLC.
 *         All rights reserved.
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
#include <cstring>

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
void tstdbcsettersandgetters( UnitTest & unitTest, int argc, char *argv[] )
{
    std::cout << "Testing Design-by-Contract setters and getters "
              << "for the UnitTest class..." << std::endl;
    
    // Silent version.
    ostringstream messages;
    
    // DBC = 0 (all off)
    {
        ScalarUnitTest foo( argc, argv, release, messages );
        foo.dbcRequire( false );
        foo.dbcCheck(   false );
        foo.dbcEnsure(  false );
        if( foo.dbcRequire() ) ITFAILS;
        if( foo.dbcCheck() )   ITFAILS;
        if( foo.dbcEnsure() )  ITFAILS;
        if( foo.dbcOn() )      ITFAILS;
    }
    // DBC = 1 (Require only)
    {
        ScalarUnitTest foo( argc, argv, release, messages );
        foo.dbcRequire( true );
        foo.dbcCheck(   false );
        foo.dbcEnsure(  false );
        if( ! foo.dbcRequire() ) ITFAILS;
        if( foo.dbcCheck() )     ITFAILS;
        if( foo.dbcEnsure() )    ITFAILS;
        if( ! foo.dbcOn() )      ITFAILS;
    }
    // DBC = 2 (Check only)
    {
        ScalarUnitTest foo( argc, argv, release, messages );
        foo.dbcRequire( false );
        foo.dbcCheck(   true );
        foo.dbcEnsure(  false );
        if( foo.dbcRequire() ) ITFAILS;
        if( ! foo.dbcCheck() ) ITFAILS;
        if( foo.dbcEnsure() )  ITFAILS;
        if( ! foo.dbcOn() )    ITFAILS;
    }
    // DBC = 4 (Ensure only)
    {
        ScalarUnitTest foo( argc, argv, release, messages );
        foo.dbcRequire( false );
        foo.dbcCheck(   false );
        foo.dbcEnsure(  true );
        if( foo.dbcRequire() )  ITFAILS;
        if( foo.dbcCheck() )    ITFAILS;
        if( ! foo.dbcEnsure() ) ITFAILS;
        if( ! foo.dbcOn() )     ITFAILS;
    }

    if( unitTest.numPasses > 0 && unitTest.numFails == 0 )
        PASSMSG( "UnitTest Design-by-Contract setters and getters are working.");
    
    return;
}

//---------------------------------------------------------------------------------------//
void tstVersion(UnitTest &ut, char *test)
{
    // Check version construction

    char version[strlen("--version")+1];
    strcpy(version, "--version");
    char *ptr1 = version;
    char *pptr[3];
    pptr[0] = test;
    pptr[2] = ptr1;
    char argument[2];
    argument[0] = 'a';
    char *ptr2 = argument;
    pptr[1] = ptr2;
    
    int argc = 3;
    char **argv = pptr;
    try
    {
        ScalarUnitTest(argc, argv, release);
        ut.failure("version construction NOT correct");
    }
    catch (assertion &err)
    {
        if (!strcmp(err.what(), "Success"))
        {
            ut.passes("version construction correct");
        }
        else
        {
            ut.failure("version construction NOT correct");
        }
    }
    catch (...)
    {
        ut.failure("version construction NOT correct");
    }
}

//---------------------------------------------------------------------------//

int main( int argc, char *argv[] )
{
    try
    {
        // Test ctor for ScalarUnitTest (also tests UnitTest ctor and member
        // function setTestName).
        ScalarUnitTest ut( argc, argv, release );

        // Also try to print the copyright and author list
        std::cout << copyright() << std::endl;

        tstOne(ut);

        // Silent version.
        ostringstream messages;
        ScalarUnitTest sut( argc, argv, release, messages );
        tstTwo(sut);

        tstTwoCheck( ut, messages );

        tstGetWordCountFile( ut );
        
        tstdbcsettersandgetters( ut, argc, argv );

        tstVersion(ut, argv[0]);
    }
    catch( rtt_dsxx::assertion &err )
    {
        string msg = err.what();
        if( msg != string( "Success" ) )
        {
            cout << "ERROR: While testing " << argv[0] << ", "
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
