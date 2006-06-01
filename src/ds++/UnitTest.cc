//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/UnitTest.cc
 * \author Kelly Thompson
 * \date   Thu May 18 15:46:19 2006
 * \brief  Implementation file for UnitTest.
 * \note   Copyright 2006 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>

#include "UnitTest.hh"

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for UnitTest object.
 * \param argc The number of command line arguments provided to main.
 * \param argv A list of command line arguments.
 * \param release_ A function pointer to the local package's release()
 * function.
 * \param out_ A user selectable output stream.  By default this is
 * std::cout. 
 *
 * This constructor automatically parses the command line to setup the name of
 * the unit test (used when generating status reports).  The object produced
 * by this constructor will respond to the command line argument "--version."
 */
UnitTest::UnitTest( int &argc, char **&argv, string_fp_void release_,
                    std::ostream & out_ )
    : testName( setTestName( std::string(argv[0])) ),
      release(   release_ ),
      numPasses( 0 ),
      numFails(  0 ),
      out( out_ )
{
    return;
}

//---------------------------------------------------------------------------//
//! Build the final message that will be desplayed when UnitTest is destroyed. 
std::string UnitTest::resultMessage()
{
    std::ostringstream msg;
    msg << "\n*********************************************\n";
    if( UnitTest::numPasses > 0 && UnitTest::numFails == 0 ) 
        msg << "**** " << testName << " Test: PASSED.\n";
    else
        msg << "**** " << testName << " Test: FAILED.\n";
    msg << "*********************************************\n";
    
    return msg.str();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Helper function to strip path information from the filename.
 * \param fqName A fully qualified filename (/path/to/the/unit/test)
 *
 * This function expects a fully qualfied name of a unit test (e.g.:
 * argv[0]).  It strips off the path and returns the name of the unit test. 
 */
std::string UnitTest::setTestName( std::string const fqName )
{
    using std::string;
    string::size_type idx=fqName.rfind('/');
    if( idx == string::npos )
        return fqName;
    string shortName = fqName.substr(idx+1);    
    return shortName;
}

//---------------------------------------------------------------------------//
/*!\brief Increment the failure count and print a message with the source line
 * number.
 * \param line The line number of the source code where the failure was
 * ecnountered. 
 */
bool UnitTest::failure(int line)
{
    out << "Test: failed on line " << line << std::endl;
    UnitTest::numFails++;
    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Increment the failure count and print a message with the source line
 * number. 
 * \param line The line number of the source code where fail was called from.
 * \param file The name of the file where the failure occured.
 */
bool UnitTest::failure(int line, char *file)
{
    out << "Test: failed on line " << line << " in " << file
	      << std::endl;
    UnitTest::numFails++;
    return false;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Increment the pass count and print a message that a test passed.
 * \param passmsg The message to be printed to the iostream \c UnitTest::out.
 */
bool UnitTest::passes(const std::string &passmsg)
{
    out << "Test: passed" << std::endl;
    out << "     " << passmsg << std::endl;
    UnitTest::numPasses++;
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Increment the failure count and print a message that a test failed.
 * \param failmsg The message to be printed to the iostream \c UnitTest::out.
 */
bool UnitTest::failure(const std::string &failmsg)
{
    out << "Test: failed" << std::endl;
    out << "     " << failmsg << std::endl;
    UnitTest::numFails++;
    return false;
}

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
//                 end of UnitTest.cc
//---------------------------------------------------------------------------//
