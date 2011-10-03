//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstGandolf_Interpreter.cc
 * \author Kelly Thompson
 * \date   Mon Oct 3 13:16 2011
 * \brief  Execute the binary Gandolf_Interpreter by redirecting the
 *         contents of GandolfInterpreter.stdin as stdin.
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "config.h"
#include "cdi_gandolf_test.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cerrno>
#include <cstring>

using namespace std;


//---------------------------------------------------------------------------//
// In this unit test we need to check binary utility Gandolf_Interpreter's
// ability to accept a collection of commands from standard input. 
// ---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void runtest()
{
    // String to hold command that will start the test.  For example:
    // "../bin/Gandolf_Interpreter Al_BeCu.ipcress < GandolfInterpreter.stdin"
    std::ostringstream cmd;
    cmd << GANDOLF_INTERPRETER_BIN_DIR << "/Gandolf_Interpreter"
        << " " << GANDOLF_INTERPRETER_BUILD_DIR << "/Al_BeCu.ipcress"
        << " < " << GANDOLF_INTERPRETER_SOURCE_DIR << "/GandolfInterpreter.stdin" 
        << " > " << GANDOLF_INTERPRETER_BUILD_DIR << "/tstGandolfInterpreter.out";
    std::string consoleCommand( cmd.str() );

    // return code from the system command
    int errorLevel(-1);
    
    // run the test.
    errno = 0;
    std::cout << consoleCommand.c_str() << std::endl;
    errorLevel = system( consoleCommand.c_str() );
    
    // check the errorLevel
    std::ostringstream msg;
    // On Linux, the child process information is sometimes unavailable even
    // for a correct run.
    if( errorLevel == 0 || errno == ECHILD)
    {
        msg << "Successful execution of Gandolf_Interpreter:"
            << "\n\t Standard input from: Gandolf_Interpreter.stdin\n";
        PASSMSG( msg.str() );
    }
    else
    {
        msg << "Unsuccessful execution of Gandolf_Interpreter:"
            << "\n\t Standard input from: Gandolf_Interpreter.stdin\n"
            << "\t errorLevel = " << errorLevel << endl
            << "\t errno = " << errno << ", " << strerror(errno) << endl;
        FAILMSG( msg.str() );   
    }
    
    // This setup provides the possibility to parse output from each
    // run.  This potential feature is not currently implemented.            
    
    return;
}

//---------------------------------------------------------------------------//
void check_output()
{
    // file contents
    std::ostringstream data; // consider using UnitTest::get_word_count(...)
    std::vector<std::string> dataByLine;

    // open the file and extract contents
    {
        std::string filename( std::string( GANDOLF_INTERPRETER_BUILD_DIR )
                              + string("/tstGandolfInterpreter.out") );
        std::ifstream infile;
        infile.open( filename.c_str() );
        Insist( infile, std::string("Cannot open specified file = \"")
                + filename + std::string("\".") );
        
        // read and store the text file contents
        std::string line;
        if( infile.is_open() )
            while( infile.good() )
            {
                getline(infile,line);
                data << line << std::endl;
                dataByLine.push_back( line );
            }
        
        infile.close();
    }

    // Spot check the output.
    {
        std::cout << "Checking the generated output file..." << std::endl;
        
        if( dataByLine[1] != std::string("This opacity file has 2 materials:") )
            ITFAILS;
        if( dataByLine[2] != std::string("Material 1 has ID number 10001") )
            ITFAILS;
        if( dataByLine[24] != std::string("Frequency grid") )
            ITFAILS;
        if( dataByLine[25] != std::string("1	1.0000000e-05") )
            ITFAILS;
        if( dataByLine[68] != std::string("material 0 Id(10001) at density 1.0000000e-01, temperature 1.0000000e+00 is: ") )
            ITFAILS;
        if( dataByLine[69] != std::string("Index 	 Group Center 		 Opacity") )
            ITFAILS;
        if( dataByLine[70] != std::string("1	 8.9050000e-03   	 1.0000000e+10") )
            ITFAILS;
    }
    
    return;
}

//---------------------------------------------------------------------------//

int main(int /*argc*/, char* /*argv*/ [])
{
    try
    {
        runtest();
        if( rtt_cdi_gandolf_test::passed )
            check_output();
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing tstGandolf_Interpreter.cc, " 
             << err.what() << endl;
        return 1;
    }
    catch( ... )
    {
        cout << "ERROR: While testing tstGandolf_Interpreter.cc, " 
             << "An unknown exception was thrown. "<< endl;
        return 1;
    }

    // status of test
    cout <<     "\n*********************************************" ;
    if( rtt_cdi_gandolf_test::passed )
        cout << "\n**** tstGandolf_Interpreter.cc Test: PASSED";
    else
        cout << "\n**** tstGandolf_Interpreter.cc Test: FAILED";
    cout <<     "\n*********************************************"
         << "\n\nDone testing tstGandolf_Interpreter.cc." << endl;

    return 0;
}   

//---------------------------------------------------------------------------//
//      end of tstGandolf_Interpreter.cc
//---------------------------------------------------------------------------//
