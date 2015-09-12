//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/ScalarUnitTest.cc
 * \author Kelly Thompson
 * \date   Thu May 18 17:08:54 2006
 * \brief  Provide services for scalar unit tests.
 * \note   Copyright (C) 2006-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ScalarUnitTest.hh"
#include "Assert.hh"
#include <iostream>
#include <sstream>
#include "XGetopt.hh"

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

    Require( argc > 0 );
    Require( release != NULL );

    // header
    out << "\n============================================="
        << "\n=== Scalar Unit Test: " << testName
        << "\n=============================================\n" << endl;

    // version tag
    out << testName << ": version " << release() << "\n" << endl;

    // exit if command line contains "--version"
    rtt_dsxx::optind=1; // resets global counter (see XGetopt.cc)

    std::map< std::string, char> long_option;
    long_option["version"] = 'v';

    int c(0);
    while(( c = rtt_dsxx::xgetopt (argc, argv, (char*)"v:", long_option)) != -1)
    {
        switch (c)
        {
            case 'v': // --version
            {
                throw rtt_dsxx::assertion( string( "Success" ) );
                return;
            }

            default:
                break; // nothing to do.
        }
    }

    Ensure( numPasses == 0 );
    Ensure( numFails  == 0 );
    Ensure( testName.length() > 0 );

    return;
}

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of ScalarUnitTest.cc
//---------------------------------------------------------------------------//
