//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/ApplicationUnitTest.cc
 * \author Kelly Thompson
 * \date   Thu Jun  1 17:15:05 2006
 * \brief  Implementation file for encapsulation of Draco application unit tests.
 * \note   Copyright (C) 2006-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ApplicationUnitTest.hh"
#include "c4/config.h"
#include "ds++/path.hh"
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <fstream>

namespace rtt_c4
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for ApplicationUnitTest
 * \arg argc The number of command line arguments
 * \arg argv A list of strings containg the command line arguments
 * \arg release_ A function pointer to this package's release function.
 * \arg applicationName A string that provides the relative path to the
 * application to be tested.
 * \arg listOfArgs_, A list of command line arguments to be provided upon
 * execution of the application.  (This might be an input deck name.)
 * \arg out_ A user specified iostream that defaults to std::cout.
 * \exception rtt_dsxx::assertion An exception with the message "Success" will
 * be thrown if \c --version is found in the argument list.  
 *
 * The constructor initializes the base class UnitTest by setting numPasses
 * and numFails to zero.  It also prints a message that declares this to be a
 * scalar unit test and provides the unit test name.
 */
ApplicationUnitTest::ApplicationUnitTest(
    int    & argc,
    char **& argv,
    string_fp_void                   release_,
    std::string              const   applicationName_,
    std::list< std::string > const & listOfArgs_,
    std::ostream & out_ )
    : UnitTest( argc, argv, release_, out_ ),
      applicationName( getFilenameComponent(applicationName_, rtt_dsxx::FC_NAME ) ),
      applicationPath( getFilenameComponent(applicationName_, rtt_dsxx::FC_PATH ) ),
      numProcs( getNumProcs( argc, argv ) ),
      mpiCommand( constructMpiCommand( numProcs ) ),
      logExtension( buildLogExtension( numProcs ) ),
      listOfArgs( listOfArgs_ ),
      logFile(),
      reportTimings(false)
{
    using std::string;

    // We need the name of the unit test, the command line argument "--np" and
    // and integer that represents the number of processors to use.
    Insist( argc > 2, "Application Tests require the argument --np [<NNN>|scalar]" );
    Require( release != NULL );
    Require( applicationName.length() > 0 );
    Require( numProcs == string("scalar") ||
             numProcs == std::string("serial") ||
             std::atoi( numProcs.c_str() ) > 0 );
    Require( mpiCommand.length() > 0 );
    Require( logExtension.length() > 0 );
    Require( testName.length() > 0 );
    
    // header
    
    out << "\n============================================="
        << "\n=== Application Unit Test: " << testName
        << "\n=== Number of Processors: "  << numProcs
        << "\n=============================================\n"
        << std::endl;
    
    // version tag
    
    out << testName << ": version " << release() << "\n" << std::endl;
    
    // exit if command line contains "--version"
    
    for( int arg = 1; arg < argc; arg++ )
    {
        if( string( argv[arg] ) == "--version" )
            throw rtt_dsxx::assertion( string( "Success" ) );

        if (string(argv[arg]) == "--timings")
            reportTimings = true;
    }

    Ensure( numPasses == 0 );
    Ensure( numFails  == 0 );
     
    return;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Print a summary of the pass/fail status of ApplicationUnitTest.
 */
void ApplicationUnitTest::status()
{
    out << "\nDone testing " << testName << "." << std::endl;
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine the number of processors to use.
 * \param argc Must be at least 3.
 * \param argv Must contain "--np [<NNN>|scalar]"
 * \retval numProcs A string holding an integral number or the name "scalar"
 */
std::string ApplicationUnitTest::getNumProcs( int & argc, char **&argv )
{
   Require(argc > 2);
   std::string np;
   // command line arguments
   for( int arg = 1; arg < argc; arg++ )
      if( std::string( argv[arg] ) == "--np" &&
          arg + 1 < argc )
          np = argv[++arg];
   Ensure( np == std::string("scalar") ||
           np == std::string("serial") ||
           std::atoi( np.c_str() ) > 0 );
    return np;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Obtain a filename extension.
 * \param numProcs The number of processor to be used for this test.
 *
 * Determine an filename extension to be used when logging the executing of
 * each input deck with the specified binary.  The actual extension takes the
 * form "-<np>.log", where <np> is the number of processors or the word
 * "scalar."
 *
 * This function provides the logic needed for future adaptation of user
 * selectable extension names.
 */
std::string ApplicationUnitTest::buildLogExtension( std::string const & numProcs )
{
    std::ostringstream le;
    Require( numProcs == std::string("scalar") ||
             numProcs == std::string("serial") ||
             std::atoi( numProcs.c_str() ) > 0 );
    le << "-";
    if( numProcs == "scalar" || numProcs == "serial" )
        le << "scalar";
    else
    {
        { // Extra check
            int np;
            std::stringstream ss(numProcs);
            ss >> np;
            Ensure( np > 0 );
        }
        le << numProcs;
    }
    le << ".out";
    Ensure( le.str().length() > 1 );
    return le.str(); 
}

//---------------------------------------------------------------------------//
/*!
 * \brief Determine the UNIX command that will be used to run the specified
 * binary.  In particular, determine if an MPI directive is required (mpirun
 * or prun) or if we are running in scalar mode.
 * \param numProcs The number of processors for this test: "scalar" or unsigned>0.
 */
std::string ApplicationUnitTest::constructMpiCommand(
    std::string const & numProcs )
{
    Require( numProcs == std::string("scalar") ||
             numProcs == std::string("serial") ||
             std::atoi( numProcs.c_str() ) > 0 );

#if defined( draco_isWin )
    { // The binary should exist.  Windows does not provide an execute bit.  
         std::string exeExists( applicationPath + applicationName + ".exe" );
         Require( std::ifstream( exeExists.c_str() ) );
    }
#else             
    { // The binary should exist and marked by the filesystem as executable.  
        
        std::string exeExistsAndExecutable("test -x " + applicationPath
                                           + applicationName );
        Require( std::system( exeExistsAndExecutable.c_str() ) == 0 );
    }
#endif

    // cmd will contain the UNIX command that will be used to execute the
    // specified program.  Possible formats include:
    //     ../bin/program
    //     mpirun -np 3 ../bin/program
    //     prun -n 16 ../bin/program
    std::ostringstream cmd;
    if( numProcs == "scalar" )
        cmd << applicationPath + applicationName << " "; // "../bin/serrano ";
    else
    {
        std::ostringstream testUname;
        
        // Determine system type:
        // On Linux use mpirun. On OSF1 use prun.
        // This information is set in config.h and in ApplicationUnitTest.hh.
        cmd << C4_MPICMD ;
            
        // relative path to the binary.
        cmd << numProcs << " " << applicationPath + applicationName
            + rtt_dsxx::exeExtension;
    }

    Ensure( cmd.str().length() > 0 );
    return cmd.str();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Add another command line argument to the list.
 * \param appArg A string spcifying a local input deck.
 */
void ApplicationUnitTest::addCommandLineArgument( std::string const & appArg )
{
    Require( appArg.length() > 0 );
    listOfArgs.push_back( appArg );
    return;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Run the application test using the arguments provided
 * 
 */
bool ApplicationUnitTest::runTest( std::string const & appArg ) 
{
    std::ostringstream msg;
    int errorLevel(0);

    std::string seperator("");
    if( appArg.length() > 0 ) seperator = std::string("_");
    logFile = testPath + applicationName + seperator
              + appArg + logExtension;
    
    std::ostringstream unixCommand;
    unixCommand << mpiCommand << " " << appArg << " > " << logFile;
    std::cout << "\nExecuting command from the shell: \n\t\""
              << unixCommand.str().c_str() << "\"\n" << std::endl;
    if (reportTimings)
    {
        problemTimer.reset();
        problemTimer.start();
    }
    errorLevel = std::system( unixCommand.str().c_str() );
    if (reportTimings)
        problemTimer.stop();

    bool result(false);
    if( errorLevel == 0 )
    {
        msg << "Test: passed\n\tSuccessful ";
        this->numPasses++;
        result = true;
    }
    else
    {
        msg << "Test: failed\n\tUnsuccessful ";
        this->numFails++;
        result = false;
    }
    msg << "execution of " << testPath + applicationName << " :"
        << "\n\tArguments           : " << appArg
        << "\n\tOutput Logfile      : " << logFile << std::endl;
    if (reportTimings)
        msg << "\tWall time           : " << problemTimer.sum_wall_clock()
            << std::endl;
    out << msg.str() << std::endl;
    return result;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Run all of the tests using all of the provided argument lists.
 * \bug Need to remove spaces and other invalid characters from log filename.
 */
void ApplicationUnitTest::runTests() 
{
    // Loop over list of input files.  Run specified binary with each, saving
    // the output to a log file.
    
    for( std::list< std::string >::const_iterator it_arg=listOfArgs.begin();
         it_arg != listOfArgs.end();
         ++it_arg )
    {
        runTest( *it_arg );
    }
    return;
}


} // end namespace rtt_c4

//---------------------------------------------------------------------------//
//                 end of ApplicationUnitTest.cc
//---------------------------------------------------------------------------//
