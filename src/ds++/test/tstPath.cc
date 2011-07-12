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

#include "../Assert.hh"
#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include "../path.hh"
#include <cstdlib> // system()

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void test_currentPath( ScalarUnitTest & ut )
{
    // currentPath
    string const cp = currentPath();

    // if we got here, currentPath didn't throw.
    
    // This path should the relative path between the top level build
    // directory and the ds++ build directory.
    string const expPartialPath( "draco/src/ds++" );

    // look for the expected partial path
    string::size_type idx = cp.find( expPartialPath );
    if( idx != string::npos )
        ut.passes( string("Found expected partial path. Path = ") + cp);
    else
        ut.failure("Did not find expected partial path. Expected path = "
                   + cp);
    
    return;
}

//---------------------------------------------------------------------------//
void test_getFilenameComponent( ScalarUnitTest & ut, string const & fqp )
{
    bool usesUnixDirSep=true;
    
    // test the FC_PATH mode
    // ------------------------------------------------------------
    
    string mypath = getFilenameComponent( fqp, rtt_dsxx::FC_PATH );
    
    string::size_type idx = mypath.find( string("test") + rtt_dsxx::UnixDirSep );
    // 2nd chance with alternate dirsep char
    if( idx == string::npos )
    {
        idx = mypath.find( string("test") + rtt_dsxx::WinDirSep );
        usesUnixDirSep = false;
    }   
    if( idx != string::npos )
        ut.passes( string("Found expected partial path. Path = ") + mypath );
    else
        ut.failure("Did not find expected partial path. Expected path = "
                   + mypath );

    // value if not found

    string mypath2 = getFilenameComponent( string("foobar"),
                                           rtt_dsxx::FC_PATH );
    if( mypath2 == string("./") )
        ut.passes("FC_PATH: name w/o path successfully returned ./");
    else
        ut.failure("FC_PATH: name w/o path returned incorrect value.");
    
    
    // test the FC_NAME mode
    // ------------------------------------------------------------
    
    string myname = getFilenameComponent( fqp, rtt_dsxx::FC_NAME );
    
    idx = myname.find( string("tstPath") );
    if( idx != string::npos )
        ut.passes( string("Found expected filename. myname = ") + myname );
    else
        ut.failure("Did not find expected filename. Expected filename = "
                   + myname );
    
    if( mypath+myname == fqp )
        ut.passes( string("Successfully divided fqp into path+name = ")
            + mypath+myname );
    else
        ut.failure( "mypath+myname != fqp" );


    // value if not found

    string myname2 = getFilenameComponent( string("foobar"),
                                           rtt_dsxx::FC_NAME );
    if( myname2 == string("foobar") )
        ut.passes("name w/o path successfully returned");
    else
        ut.failure("name w/o path returned incorrect value.");
        

    // test the FC_REALPATH
    // ------------------------------------------------------------

    string realpath = getFilenameComponent( fqp, rtt_dsxx::FC_REALPATH );
    
#if defined( draco_isWin )
    { // The binary should exist.  Windows does not provide an execute bit.  
         std::string exeExists( realpath + ".exe" );
         if( std::ifstream( exeExists.c_str() ) )
            ut.passes( "FC_REALPATH points to a valid executable." );
        else
            ut.failure( "FC_REALPATH is invalid or not executable." );
    }
#else             
    { // The binary should exist and marked by the filesystem as executable.  

        if( usesUnixDirSep )
        {
            if( realpath.size() > 0 ) ut.passes(  "FC_REALPATH has length > 0.");
            else                      ut.failure( "FC_REALPATH has length <= 0.");

            std::string exeExistsAndExecutable("test -x " + realpath );
            if( std::system( exeExistsAndExecutable.c_str() ) == 0 )
                ut.passes(
                    string("FC_REALPATH points to a valid executable. Path = ")
                    + realpath );
            else
                ut.failure(
                    string( "FC_REALPATH is invalid or not executable. Path = ")
                    + realpath );
        }
        else
        {
            if( realpath.size() == 0 ) ut.passes(  "FC_REALPATH has length <= 0.");
            else                       ut.failure( "FC_REALPATH has length > 0.");
            
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
        ut.passes( "FC_ABSOLUTE throws." );
    }
    if( ! caught ) ut.failure( "FC_ABSOLUTE failed to throw." );

    caught=false;
    try
    {
        string extension = getFilenameComponent( fqp, rtt_dsxx::FC_EXT );
    }
    catch( ... )
    {
        caught=true;
        ut.passes( "FC_EXT throws." );
    }
    if( ! caught ) ut.failure( "FC_EXT failed to throw." );
    
    caught=false;
    try
    {
        string extension = getFilenameComponent( fqp, rtt_dsxx::FC_NAME_WE );
    }
    catch( ... )
    {
        caught=true;
        ut.passes( "FC_NAME_WE throws." );
    }
    if( ! caught ) ut.failure( "FC_NAME_WE failed to throw." );

    // FC_LASTVALUE should always throw.
    caught=false;
    try
    {
        string foo = getFilenameComponent( fqp, rtt_dsxx::FC_LASTVALUE );
    }
    catch( ... )
    {
        caught=true;
        ut.passes( "FC_LASTVALUE throws." );
    }
    if( ! caught ) ut.failure( "FC_LASTVALUE failed to throw." );
    
    return;
}


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_dsxx::ScalarUnitTest ut(argc, argv, release);
    try
    {
        test_currentPath(ut);
        test_getFilenameComponent(
            ut,
            string("test") + rtt_dsxx::UnixDirSep + string("tstPath"));
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
