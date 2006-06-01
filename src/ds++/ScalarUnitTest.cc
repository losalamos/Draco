//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/ScalarUnitTest.cc
 * \author Kelly Thompson
 * \date   Thu May 18 17:08:54 2006
 * \brief  Provide services for scalar unit tests.
 * \note   Copyright 2006 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>

#include "ScalarUnitTest.hh"
#include "Assert.hh"

namespace rtt_dsxx
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for ScalarUnitTest
 * \arg argc The number of command line arguments
 * \arg argv A list of strings containg the command line arguments
 * \arg release_ A function pointer to this package's release function.
 * \arg out_ A user specified iostream that defaults to std::cout.
 * \exception rtt_dsxx::assertion An exception with the message "Success" will
 * be thrown if \c --version is found in the argument list.  
 *
 * The constructor initializes the base class UnitTest by setting numPasses
 * and numFails to zero.  It also prints a message that declares this to be a
 * scalar unit test and provides the unit test name.
 */
ScalarUnitTest::ScalarUnitTest( int &argc, char **&argv,
                                string_fp_void release_,
                                std::ostream & out_ )
    : UnitTest( argc, argv, release_, out_ )
{
    using std::endl;
    using std::string;

    // header
    
    out << "\n============================================="
         << "\n=== Scalar Unit Test: " << testName
         << "\n=============================================\n" << endl;

    // version tag

    out << testName << ": version " << release() << "\n" << endl;

    // exit if command line contains "--version"
    
     for( int arg = 1; arg < argc; arg++ )
         if( string( argv[arg] ) == "--version" )
             throw rtt_dsxx::assertion( string( "Success" ) );
     
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor
 *
 * The destructor will provides a final status report during destruction.
 */
ScalarUnitTest::~ScalarUnitTest()
{
    out << resultMessage() << std::endl;
    return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Print a summary of the pass/fail status of UnitTest.
 */
void ScalarUnitTest::status()
{
    out << "\nDone testing " << testName  << "." << std::endl;
    return;
}

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
//                 end of ScalarUnitTest.cc
//---------------------------------------------------------------------------//
