//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/driver4tstConsole_Token_Stream.cc
 * \author Kelly Thompson
 * \date   Wed Oct 19 14:42 2005
 * \brief  Execute the binary tstConsole_Token_Stream by redirecting the
 * contents of console_test.inp as stdin.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>

#include "c4/global.hh"
#include "c4/SpinLock.hh"
#include "../Release.hh"
#include "parser_test.hh"

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
    
    std::ostringstream unixCommand;
    
    // Determine system type:
    
    const bool isLinux( system("test `uname` = Linux") == 0 );
    const bool isOSF1( system("test `uname` = OSF1") == 0 );
    
    // On Linux use mpirun. On OSF1 use prun.
    
    if( isLinux )
    {
        std::cout << "\nSystem type is Linux" << std::endl;
        unixCommand << "mpirun -np ";
    }
    else if( isOSF1 )
    {
        std::cout << "\nSystem type is OSF1" << std::endl;
        unixCommand << "prun -n ";
    }
    
    // relative path to the executable
    
    unixCommand << rtt_c4::nodes()
                << " ./tstConsole_Token_Stream < console_test.inp ";
    
    // return code from the system command
    int errorLevel(-1);
    
    // run the test.
    errorLevel = system( unixCommand.str().c_str() );
    
    // check the errorLevel
    std::ostringstream msg;
    if( errorLevel == 0 )
    {
        msg << "Successful execution of tstConsole_Token_Stream:"
            << "\n\t Number of processors: " << rtt_c4::nodes()
            << "\n\t Input deck          : console_test.inp\n";
        PASSMSG( msg.str() );
    }
    else
    {
        msg << "Unsuccessful execution of tstConsole_Token_Stream:"
            << "\n\t Number of processors: " << rtt_c4::nodes()
            << "\n\t Input deck          : console_test.inp\n";
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
        if (std::string(argv[arg]) == "--version")
        {
            rtt_c4::finalize();
            return 0;
        }

    try
    {
        runtest();
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing driver4tstConsole_Token_Stream.cc, " 
                  << err.what() << std::endl;
        rtt_c4::finalize();
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing driver4tstConsole_Token_Stream.cc, " 
                  << "An unknown exception was thrown on processor "
                  << rtt_c4::node() << std::endl;
        rtt_c4::finalize();
        return 1;
    }

    {
        rtt_c4::HTSyncSpinLock slock;

        // status of test
        std::cout << std::endl;
        std::cout <<     "*********************************************" 
                  << std::endl;
        if (rtt_parser_test::passed) 
        {
            std::cout << "**** driver4tstConsole_Token_Stream.cc Test: PASSED on " 
                      << rtt_c4::node() 
                      << std::endl;
        }
        std::cout <<     "*********************************************" 
                  << std::endl;
        std::cout << std::endl;
    }
    
    rtt_c4::global_barrier();

    std::cout << "Done testing driver4tstConsole_Token_Stream.cc on " << rtt_c4::node() 
              << std::endl;
    
    rtt_c4::finalize();

    return 0;
}   

//---------------------------------------------------------------------------//
//      end of driver4tstConsole_Token_Stream.cc
//---------------------------------------------------------------------------//
