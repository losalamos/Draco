//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstPath.cc
 * \author Kelly Thompson
 * \date   Tue Jul 12 16:00:59 2011
 * \brief  Test functions found in ds++/path.hh and path.cc
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/config.h"
#include "../Assert.hh"
#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include "../path.hh"
#include "../SystemCall.hh"
#include <cstdlib> // system()
#include <fstream>

using namespace std;
using namespace rtt_dsxx;

#define PASSMSG(m) ut.passes(m)
#define FAILMSG(m) ut.failure(m)
#define ITFAILS    ut.failure( __LINE__, __FILE__ )

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void test_currentPath( ScalarUnitTest & ut )
{
    cout << "\nTesting currentPath() function ... \n" << endl;
    
    // currentPath
    string const cp = draco_getcwd();

    // if we got here, currentPath didn't throw.
    
    // Note, we have no idea where the this test was run from so we can say
    // nothing about what the path string should contain.
    
    if( fileExists( cp ) )
        PASSMSG( string("Retrieved current path exists. cp = ") + cp);
    else
        FAILMSG("Retrieved current path does not exist. cp = " + cp);

    // Test behavior of fileExist when file does not exist.
    string const fileDoesNotExist( "/bounty/hunter/boba_fett" );
        
    if( ! fileExists( fileDoesNotExist ) )
        PASSMSG( string("fileExist() correctly returned false. ")) ;
    else
        FAILMSG( string("fileExist() incorrectly returned true. ")) ;
    
    return;
}

//---------------------------------------------------------------------------//
void test_getFilenameComponent( ScalarUnitTest & ut, string const & fqp )
{
    cout << "\nTesting getFilenameComponent() function with fqp = "
         << fqp << " ...\n" << endl;

    bool usesUnixDirSep=true;
    
    // test the FC_PATH mode
    // ------------------------------------------------------------

    // 4 possible cases: ./tstPath, .../test/tstPath,
    // tstPath.exe or test/tstPath.exe

    // Does the provided path use unix or windows directory separator?
    string::size_type idx = fqp.find( rtt_dsxx::UnixDirSep );
    if( idx == string::npos )
        usesUnixDirSep = false;

    // retrieve the path w/o the filename.
    string mypath = getFilenameComponent( fqp, rtt_dsxx::FC_PATH );

    if( usesUnixDirSep )
    {
        // If we are using UnixDirSep, we have 2 cases (./tstPath or
        // .../test/tstPath).  Look for the case with 'test' first:
        idx = mypath.find( string("test") + rtt_dsxx::UnixDirSep );

        // If the return string does not have 'test/' then we also need to
        // check for './' as a pass condition
        if( idx == string::npos )
            idx = mypath.find( rtt_dsxx::UnixDirSep );

        // Report pass/fail
        if( idx != string::npos )
            PASSMSG( string("Found expected partial path. Path = ") + mypath );
        else
            FAILMSG("Did not find expected partial path. Expected path = "
                       + mypath );
    }
    else
    {
        // If we are using WinDirSep, we have 2 cases (.\tstPath.exe or
        // ...\test\tstPath.exe).  Look for the case with 'test' first:
        idx = mypath.find( string("test") + rtt_dsxx::WinDirSep );

        // If the return string does not have 'test\' then we also need to
        // check for './' as a pass condition
        if( idx == string::npos )
            idx = mypath.find( rtt_dsxx::WinDirSep );

        // Report pass/fail
        if( idx != string::npos )
            PASSMSG( string("Found expected partial path. Path = ") + mypath );
        else
            FAILMSG("Did not find expected partial path. Expected path = "
                       + mypath );
    }
    
    // value if not found

    string mypath2 = getFilenameComponent( string("foobar"),
                                           rtt_dsxx::FC_PATH );
    string expected = string(".")+string(1,rtt_dsxx::dirSep);
    if( mypath2 == expected )
        PASSMSG( string("FC_PATH: name w/o path successfully returned ")
                   +expected);
    else
        FAILMSG("FC_PATH: name w/o path returned incorrect value.");
    
    
    // test the FC_NAME mode
    // ------------------------------------------------------------
    
    string myname = getFilenameComponent( fqp, rtt_dsxx::FC_NAME );
    
    idx = myname.find( string("tstPath") );
    if( idx != string::npos )
        PASSMSG( string("Found expected filename. myname = ") + myname );
    else
        FAILMSG("Did not find expected filename. Expected filename = "
                   + myname );
        
    if( usesUnixDirSep )
    {
        if( mypath+myname == fqp )
            PASSMSG( string("Successfully divided fqp into path+name = ")
                        + mypath+myname );
        else
            FAILMSG( string("mypath+myname != fqp") +
                        string("\n\tmypath = ") + mypath +
                        string("\n\tmyname = ") + myname +
                        string("\n\tfqp    = ") + fqp
                      );
    } 

    // value if not found

    string myname2 = getFilenameComponent( string("foobar"),
                                           rtt_dsxx::FC_NAME );
    if( myname2 == string("foobar") )
        PASSMSG("name w/o path successfully returned");
    else
        FAILMSG("name w/o path returned incorrect value.");
        

    // test the FC_REALPATH
    // ------------------------------------------------------------

    string realpath = getFilenameComponent( fqp, rtt_dsxx::FC_REALPATH );

#if defined( WIN32 )
    { // The binary should exist.  Windows does not provide an execute bit.  
         std::string exeExists( realpath + ".exe" );
         if( std::ifstream( exeExists.c_str() ) )
            PASSMSG( "FC_REALPATH points to a valid executable." );
         else
            FAILMSG( string("FC_REALPATH is invalid or not executable.") +
                     string("  exeExists = ") + exeExists );
    }
#else             
    { // The binary should exist and marked by the filesystem as executable.  

        if( usesUnixDirSep )
        {
            if( realpath.size() > 0 ) PASSMSG( "FC_REALPATH has length > 0.");
            else                      FAILMSG( "FC_REALPATH has length <= 0.");

            std::string exeExistsAndExecutable("test -x " + realpath );
            if( std::system( exeExistsAndExecutable.c_str() ) == 0 )
                PASSMSG(
                    string("FC_REALPATH points to a valid executable. Path = ")
                    + realpath );
            else
                FAILMSG(
                    string( "FC_REALPATH is invalid or not executable. Path = ")
                    + realpath );
        }
        else
        {
            if( realpath.size() == 0 ) PASSMSG( "FC_REALPATH has length <= 0.");
            else                       FAILMSG( "FC_REALPATH has length > 0.");
            
        }
    }
#endif


    // These are not implemented yet
    // ------------------------------------------------------------

    bool caught=false;
    try
    {
        string absolutepath = getFilenameComponent( fqp, rtt_dsxx::FC_ABSOLUTE );
    }
    catch( ... )
    {
        caught=true;
        PASSMSG( "FC_ABSOLUTE throws." );
    }
    if( ! caught ) FAILMSG( "FC_ABSOLUTE failed to throw." );

    caught=false;
    try
    {
        string extension = getFilenameComponent( fqp, rtt_dsxx::FC_EXT );
    }
    catch( ... )
    {
        caught=true;
        PASSMSG( "FC_EXT throws." );
    }
    if( ! caught ) FAILMSG( "FC_EXT failed to throw." );
    
    caught=false;
    try
    {
        string extension = getFilenameComponent( fqp, rtt_dsxx::FC_NAME_WE );
    }
    catch( ... )
    {
        caught=true;
        PASSMSG( "FC_NAME_WE throws." );
    }
    if( ! caught ) FAILMSG( "FC_NAME_WE failed to throw." );

    // FC_LASTVALUE should always throw.
    caught=false;
    try
    {
        string foo = getFilenameComponent( fqp, rtt_dsxx::FC_LASTVALUE );
    }
    catch( ... )
    {
        caught=true;
        PASSMSG( "FC_LASTVALUE throws." );
    }
    if( ! caught ) FAILMSG( "FC_LASTVALUE failed to throw." );
    
    return;
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_dsxx::ScalarUnitTest ut(argc, argv, release);
    try
    {
        test_currentPath(ut);

        test_getFilenameComponent(ut, string(argv[0]));
            
        test_getFilenameComponent(
            ut,
            string("test") + rtt_dsxx::WinDirSep + string("tstPath.exe"));
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing tstPath, " << err.what() << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        cout << "ERROR: While testing tstPath, " 
             << "An unknown exception was thrown." << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstPath.cc
//---------------------------------------------------------------------------//
