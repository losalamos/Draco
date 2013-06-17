//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/driver4tstConsole_Token_Stream.cc
 * \author Kelly Thompson
 * \date   Wed Oct 19 14:42 2005
 * \brief  Execute the binary tstConsole_Token_Stream by redirecting the
 *         contents of console_test.inp as stdin.
 * \note   Copyright (C) 2004-1013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/path.hh"
#include <sstream>
// #include <cstdlib>
// #include <cerrno>
#include <cstring>

using namespace std;

#define PASSMSG(a) ut.passes(a)
#define ITFAILS    ut.failure(__LINE__)
#define FAILURE    ut.failure(__LINE__, __FILE__)
#define FAILMSG(a) ut.failure(a)

//---------------------------------------------------------------------------//
// In this unit test we need to check the parser's ability to accept a
// collection of tokens from standard input.  This cannot be done with a
// single unit test.  Instead, we create the binary tstConsole_Token_Stream
// that will accept data from standard input.  In this file, we actually
// execute tstConsole_Token_Stream and pipe the data from the file
// console_test.inp as input.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void runtest( rtt_dsxx::UnitTest &ut )
{
    // Build path for the application "tstConsole_Token_Stream"

    string tctsApp = rtt_dsxx::getFilenameComponent(
        ut.getTestPath() + std::string("tstConsole_Token_Stream"),
        rtt_dsxx::FC_NATIVE) + rtt_dsxx::exeExtension;
// #ifndef _MSC_VER    
//     // Unix needs the leading dot-slash.
//     tctsApp = std::string("./") + tctsApp;
// #endif

    // Build path for the input file "console_test.inp"
    string const ctInputFile = rtt_dsxx::getFilenameComponent(
        ut.getTestPath() + std::string("console_test.inp"),
        rtt_dsxx::FC_NATIVE);

    // String to hold command that will start the test.  For example:
    // "mpirun -np 1 ./tstConsole_Token_Stream < console_test.inp"    
    std::string consoleCommand( tctsApp + " < " + ctInputFile );

    // return code from the system command
    int errorLevel(-1);
    
    // run the test.
    errno = 0;
    errorLevel = system( consoleCommand.c_str() );
    
    // check the errorLevel
    std::ostringstream msg;
    // On Linux, the child process information is sometimes unavailable even
    // for a correct run.
    if( errorLevel == 0 || errno == ECHILD)
    {
        msg << "Successful execution of tstConsole_Token_Stream:"
            << "\n\t Standard input from: console_test.inp\n";
        PASSMSG( msg.str() );
    }
    else
    {
        msg << "Unsuccessful execution of tstConsole_Token_Stream:"
            << "\n\t Standard input from: console_test.inp\n"
            << "\t errorLevel = " << errorLevel << endl
            << "\t errno = " << errno << ", " << strerror(errno) << endl;
        FAILMSG( msg.str() );   
    }
    
    // This setup provides the possibility to parse output from each
    // run.  This potential feature is not currently implemented.            
    
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_dsxx::ScalarUnitTest ut( argc, argv, rtt_dsxx::release );
    try { runtest(ut ); }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of driver4tstConsole_Token_Stream.cc
//---------------------------------------------------------------------------//
