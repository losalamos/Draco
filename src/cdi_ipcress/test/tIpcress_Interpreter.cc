//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstIpcress_Interpreter.cc
 * \author Kelly Thompson
 * \date   Mon Oct 3 13:16 2011
 * \brief  Execute the binary Ipcress_Interpreter by redirecting the
 *         contents of IpcressInterpreter.stdin as stdin.
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "config.h"
#include "cdi_ipcress_test.hh"
#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cerrno>
#include <cstring>

using namespace std;


#define PASSMSG(m) ut.passes(m)
#define FAILMSG(m) ut.failure(m)
#define ITFAILS    ut.failure( __LINE__, __FILE__ )

//---------------------------------------------------------------------------//
// In this unit test we need to check binary utility Ipcress_Interpreter's
// ability to accept a collection of commands from standard input. 
// ---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void runtest(rtt_dsxx::ScalarUnitTest &ut)
{
    // String to hold command that will start the test.  For example:
    // "../bin/Ipcress_Interpreter Al_BeCu.ipcress < IpcressInterpreter.stdin"
    ostringstream cmd;
    cmd << IPCRESS_INTERPRETER_BIN_DIR << "/Ipcress_Interpreter"
        // << " " << IPCRESS_INTERPRETER_BUILD_DIR << "/Al_BeCu.ipcress"
        << " " << "Al_BeCu.ipcress"
        << " < " << IPCRESS_INTERPRETER_SOURCE_DIR << "/IpcressInterpreter.stdin" 
        << " > " << IPCRESS_INTERPRETER_BUILD_DIR << "/tstIpcressInterpreter.out";

    cout << "Preparing to run: \n" << endl;
    string consoleCommand( cmd.str() );

    // return code from the system command
    int errorLevel(-1);
    
    // run the test.
    errno = 0;
    cout << consoleCommand.c_str() << endl;
    errorLevel = system( consoleCommand.c_str() );
    
    // check the errorLevel
    ostringstream msg;
    // On Linux, the child process information is sometimes unavailable even
    // for a correct run.
    if( errorLevel == 0 || errno == ECHILD)
    {
        msg << "Successful execution of Ipcress_Interpreter:"
            << "\n\t Standard input from: Ipcress_Interpreter.stdin\n";
        PASSMSG( msg.str() );
    }
    else
    {
        msg << "Unsuccessful execution of Ipcress_Interpreter:"
            << "\n\t Standard input from: Ipcress_Interpreter.stdin\n"
            << "\t errorLevel = " << errorLevel << endl
            << "\t errno = " << errno << ", " << strerror(errno) << endl;
        FAILMSG( msg.str() );   
    }
    
    // This setup provides the possibility to parse output from each
    // run.  This potential feature is not currently implemented.            
    
    return;
}

//---------------------------------------------------------------------------//
void check_output(rtt_dsxx::ScalarUnitTest &ut)
{
    // file contents
    ostringstream data; // consider using UnitTest::get_word_count(...)
    vector<string> dataByLine;

    // open the file and extract contents
    {
        string filename( string( IPCRESS_INTERPRETER_BUILD_DIR )
                         + string("/tstIpcressInterpreter.out") );
        ifstream infile;
        infile.open( filename.c_str() );
        Insist( infile, string("Cannot open specified file = \"")
                + filename + string("\".") );
        
        // read and store the text file contents
        string line;
        if( infile.is_open() )
            while( infile.good() )
            {
                getline(infile,line);
                data << line << endl;
                dataByLine.push_back( line );
            }
        
        infile.close();
    }

    // Spot check the output.
    {
        cout << "Checking the generated output file..." << endl;
        
        if( dataByLine[1] != string("This opacity file has 2 materials:") )
            ITFAILS;
        if( dataByLine[2] != string("Material 1 has ID number 10001") )
            ITFAILS;
        if( dataByLine[24] != string("Frequency grid") )
            ITFAILS;
        if( dataByLine[25] != string("1	1.0000000e-05") )
            ITFAILS;
        if( dataByLine[69] != string("material 0 Id(10001) at density 1.0000000e-01, temperature 1.0000000e+00 is: ") )
        {
            cout << "match failed: \n   "
                 << dataByLine[68] << " != \n   "
                 << "material 0 Id(10001) at density 1.0000000e-01, temperature 1.0000000e+00 is: " << endl;
            ITFAILS;
        }
        if( dataByLine[70] != string("Index 	 Group Center 		 Opacity") )
            ITFAILS;
        if( dataByLine[71] != string("1	 8.9050000e-03   	 1.0000000e+10") )
            ITFAILS;
    }
    
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char* argv[])
{
    rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
    try
    {
        runtest(ut);
        if( ut.numFails == 0 )
            check_output(ut);
    }
    catch ( exception const & err)
    {
        cout << "ERROR: While testing tstIpcress_Interpreter.cc, " 
             << err.what() << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        cout << "ERROR: While testing tstIpcress_Interpreter.cc, " 
             << "An unknown exception was thrown. "<< endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//      end of tstIpcress_Interpreter.cc
//---------------------------------------------------------------------------//
