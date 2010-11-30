//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/driver4tstConsole_Token_Stream.cc
 * \author Kelly Thompson
 * \date   Wed Oct 19 14:42 2005
 * \brief  Execute the binary tstConsole_Token_Stream by redirecting the
 * contents of console_test.inp as stdin.
 * \note   Copyright (C) 2004-1020 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "parser_test.hh"
#include "../Release.hh"
#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cerrno>
#include <cstring>

using namespace std;


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

void runtest()
{
    // Only initiate test from processor 0.
    if( rtt_c4::node() != 0 ) return;

    // String to hold command that will start the test.  For example:
    // "mpirun -np 1 ./tstConsole_Token_Stream < console_test.inp"
    
    std::string consoleCommand("tstConsole_Token_Stream < console_test.inp ");
#ifndef _MSC_VER    
    // Unix needs the leading dot-slash.
    consoleCommand = std::string("./") + consoleCommand;
#endif
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
    
    rtt_c4::initialize(argc, argv);

    // version tag
    if (rtt_c4::node() == 0)
        cout << argv[0] << ": version " << rtt_parser::release() << endl;

    // Optional exit
    for (int arg = 1; arg < argc; arg++)
        if (string(argv[arg]) == "--version")
        {
            rtt_c4::finalize();
            return 0;
        }

    try
    {
        runtest();
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing driver4tstConsole_Token_Stream.cc, " 
             << err.what() << endl;
        rtt_c4::finalize();
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing driver4tstConsole_Token_Stream.cc, " 
             << "An unknown exception was thrown on processor "
             << rtt_c4::node() << endl;
        rtt_c4::finalize();
        return 1;
    }

    { // status of test
        rtt_c4::HTSyncSpinLock slock;
        if (rtt_c4::node() == 0)
        {
            cout <<     "\n*********************************************" ;
            if (rtt_parser_test::passed) 
                cout << "\n**** driver4tstConsole_Token_Stream.cc Test: PASSED on " 
                << rtt_c4::node();
            cout <<     "\n*********************************************\n"
              << endl;
        }
    }
    
    rtt_c4::global_barrier();

    cout << "Done testing driver4tstConsole_Token_Stream.cc on " << rtt_c4::node() 
         << endl;
    rtt_c4::finalize();

    return 0;
}   

//---------------------------------------------------------------------------//
//      end of driver4tstConsole_Token_Stream.cc
//---------------------------------------------------------------------------//
