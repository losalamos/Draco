//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/test/tstIpcress_Interpreter.cc
 * \author Kelly Thompson
 * \date   Mon Oct 3 13:16 2011
 * \brief  Execute the binary Ipcress_Interpreter by redirecting the
 *         contents of IpcressInterpreter.stdin as stdin.
 * \note   Copyright (C) 2011-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_ipcress_test.hh"
#include "ds++/Release.hh"
#include <fstream>
#include <cstdlib>
#include <cerrno>
#include <cstring>

using namespace std;

//---------------------------------------------------------------------------//
// In this unit test we need to check binary utility Ipcress_Interpreter's
// ability to accept a collection of commands from standard input.
// ---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void runtest(rtt_dsxx::ScalarUnitTest &ut)
{
    // Setup paths:

    // The build system knows the location of Ipcress_Interpreter, it is provided
    // by CPP processor definition.
    string const ii_application( IILOC );
    string const tpath( ut.getTestPath() );       // path to test binary dir
    string const tip(   ut.getTestSourcePath() ); // path to test source dir

    // String to hold command that will start the test.  For example:
    // "../Ipcress_Interpreter Al_BeCu.ipcress < IpcressInterpreter.stdin"
    ostringstream cmd;
    cmd << ii_application << " Al_BeCu.ipcress"
        << " < " << tip   << "IpcressInterpreter.stdin"
        << " > " << tpath << "tstIpcressInterpreter.out";

    // return code from the system command
    int errorLevel(-1);

    // run the test.
    cout << "Preparing to run: \n" << cmd.str() << endl;
    errno = 0;
    errorLevel = system( cmd.str().c_str() );

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
        string filename( ut.getTestPath()
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
        // dataByLine[25] = "1\t1.0000000e-005"
        std::vector<std::string> wordsInLine
            = rtt_dsxx::UnitTest::tokenize( dataByLine[25], "\t" );
        double ddata = atof( wordsInLine[1].c_str() );
        if( ! rtt_dsxx::soft_equiv( ddata, 1.0e-5 ) )
            ITFAILS;
        // dataByLine[69] = "material 0 Id(10001) at density 1.0000000e-01, temperature 1.0000000e+00 is: "
        wordsInLine = rtt_dsxx::UnitTest::tokenize( dataByLine[69] );
        if( wordsInLine[0] != string("material") )  ITFAILS;
        if( wordsInLine[1] != string("0") )         ITFAILS;
        if( wordsInLine[2] != string("Id(10001)") ) ITFAILS;
        if( wordsInLine[3] != string("at") )        ITFAILS;
        if( wordsInLine[4] != string("density") )   ITFAILS;
        if( ! rtt_dsxx::soft_equiv( atof(
                                        wordsInLine[5].substr(0,wordsInLine[5].size()-1).c_str() ),
                                    1.0e-01 ) ) ITFAILS;
        if( wordsInLine[6] != string("temperature") )  ITFAILS;
        if( ! rtt_dsxx::soft_equiv( atof( wordsInLine[7].c_str() ), 1.0 ) )
            ITFAILS;
        if( wordsInLine[8] != string("is:") )          ITFAILS;

        if( dataByLine[70] != string("Index \t Group Center \t\t Opacity") )
            ITFAILS;
        // dataByLine[71] = "1\t8.9050000e-03\t1.0000000e+10"
        wordsInLine = rtt_dsxx::UnitTest::tokenize( dataByLine[71], "\t" );
        if( atoi( wordsInLine[0].c_str() ) != 1 ) ITFAILS;
        if( ! rtt_dsxx::soft_equiv( atof( wordsInLine[1].c_str() ), 8.9050000e-03 ) )
            ITFAILS;
        if( ! rtt_dsxx::soft_equiv( atof( wordsInLine[2].c_str() ), 1.0000000e+10 ) )
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
    UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstIpcress_Interpreter.cc
//---------------------------------------------------------------------------//
