//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/UnitTest.cc
 * \author Kelly Thompson
 * \date   Thu May 18 15:46:19 2006
 * \brief  Implementation file for UnitTest.
 * \note   Copyright (C) 2006-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "UnitTest.hh"
#include "path.hh"
#ifdef DRACO_DIAGNOSTICS_LEVEL_3
#include "fpe_trap.hh"
#endif
#include <fstream>
#include <sstream>

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
UnitTest::UnitTest( int              & /* argc */, 
                    char           **& argv,
                    string_fp_void     release_,
                    std::ostream     & out_ )
    : numPasses( 0 ),
      numFails(  0 ),
      fpe_trap_active(false),
      testName( getFilenameComponent( std::string(argv[0]), rtt_dsxx::FC_NAME )),
      testPath( getFilenameComponent( std::string(argv[0]), rtt_dsxx::FC_PATH )),
      release(  release_ ),
      out(      out_ ),
      m_dbcRequire(false),
      m_dbcCheck(false),
      m_dbcEnsure(false),
      m_dbcNothrow(false)
{
    Require( release   != NULL );
    Ensure(  numPasses == 0 );
    Ensure(  numFails  == 0 );
    Ensure(  testName.length() > 0 );
#if DBC & 1
    m_dbcRequire = true;
#endif
#if DBC & 2
    m_dbcCheck = true;
#endif
#if DBC & 4
    m_dbcEnsure = true;
#endif
#if DBC & 8
    m_dbcNothrow = true;
#endif

    // Turn on fpe_traps at level 3.
#ifdef DRACO_DIAGNOSTICS_LEVEL_3
    // if set to false, most compilers will provide a stack trace.
    // if set to true, fpe_trap forms a simple message and calls Insist.
    bool const abortWithInsist(true);
    rtt_dsxx::fpe_trap fpeTrap(abortWithInsist);
    fpe_trap_active = fpeTrap.enable();
#endif
    
    return;
}

//---------------------------------------------------------------------------//
//! Build the final message that will be desplayed when UnitTest is destroyed. 
std::string UnitTest::resultMessage() const
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
bool UnitTest::failure(int line, char const *file)
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

//---------------------------------------------------------------------------//
/*!
 * \brief Parse msg to provide a list of words and the number of occurances of each.
 */
std::map< std::string, unsigned >
UnitTest::get_word_count( std::ostringstream const & msg, bool verbose )
{
    using std::map;
    using std::string;
    using std::cout;
    using std::endl;
    
    map<string,unsigned> word_list;
    string msgbuf( msg.str() );
    string delims(" \n\t:,.;");
    
    { // Build a list of words found in msgbuf.  Count the number of
      // occurances.
        
        // Find the beginning of the first word.
        string::size_type begIdx = msgbuf.find_first_not_of(delims);
        string::size_type endIdx;
        
        // While beginning of a word found
        while( begIdx != string::npos )
        {
            // search end of actual word
            endIdx = msgbuf.find_first_of( delims, begIdx );
            if( endIdx == string::npos)
                endIdx = msgbuf.length();
            
            // the word is we found is...
            string word( msgbuf, begIdx, endIdx-begIdx );
            
            // add it to the map
            word_list[ word ]++;
            
            // search to the beginning of the next word
            begIdx = msgbuf.find_first_not_of( delims, endIdx );        
        }
    }

    if( verbose )
    {
        cout << "The messages from the message stream contained the following words/occurances." << endl;
        // print the word_list
        for( auto it = word_list.begin(); it != word_list.end(); ++it)
            cout << it->first << ": " << it->second << endl;
    }

    return word_list;
}

//---------------------------------------------------------------------------------------//
//! \brief Convert a string into a vector of words.
std::vector<std::string> UnitTest::tokenize(
    std::string const & source,
    char        const * delimiter_list,
    bool                keepEmpty)
{
    std::vector<std::string> results;

    size_t prev = 0;
    size_t next = 0;

    while ((next = source.find_first_of(delimiter_list, prev)) != std::string::npos)
    {
        if (keepEmpty || (next - prev != 0))
            results.push_back(source.substr(prev, next - prev));
        prev = next + 1;
    }

    if (prev < source.size())
        results.push_back(source.substr(prev));

    return results;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Parse text file to provide a list of words and the number of occurances of each.
 */
std::map< std::string, unsigned >
UnitTest::get_word_count( std::string const & filename, bool verbose )
{
    // open the file
    std::ifstream infile;
    infile.open( filename.c_str() );
    Insist( infile, std::string("Cannot open specified file = \"") + filename
            + std::string("\".") );

    // read and store the text file contents
    std::ostringstream data;
    std::string line;
    if( infile.is_open() )
        while( infile.good() )
        {
            getline(infile,line);
            data << line << std::endl;
        }

    infile.close();
    return UnitTest::get_word_count( data, verbose );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Returns the expected path of the input file
 *
 * IDEs often have "configuration" subdirectories like "Debug" and "Release"
 * that they place the test executables in, leaving the input files one level
 * up.  This method attempts to detect that condition in order to provide
 * the correct path, regardless of the selected build system.
 */
std::string
UnitTest::getTestInputPath() const
{
    std::string inputDir(this->getTestPath());

    // If inputDir is a relative path that points to "./" replace it with a full
    // path.
    if( inputDir == std::string("./") )
        inputDir = rtt_dsxx::draco_getcwd();

    // Next, check to see if we are working with a configuration subdirectory
    if( rtt_dsxx::getFilenameComponent( inputDir, rtt_dsxx::FC_NAME ) == "test")
    {
        // Then there is no configuration subdirectory (like "Debug") at the end
        // of the path; simply return the inputDir from the current dir
        // (This simply ensures that the slashes point the right way).
        inputDir = rtt_dsxx::getFilenameComponent(inputDir,rtt_dsxx::FC_NATIVE);
    }
    else
    {
        // The project generator does use $(configuration).  So the input
        // file is at [package]/test/file.inp, but the test binary is
        // at [package]/test/$(Configuration)/; return the former
        std::string configuration = rtt_dsxx::getFilenameComponent(
                                              inputDir, rtt_dsxx::FC_NAME );
        int pos(inputDir.find_last_of(configuration) - configuration.length());
        Check( pos > 0 );
        std::string inputDirTrunc( inputDir.substr(0,pos) + rtt_dsxx::dirSep );

        // Ensure that the slashes are correct
        inputDir = rtt_dsxx::getFilenameComponent(inputDirTrunc,
                                                  rtt_dsxx::FC_NATIVE);
    }

    return (inputDir);
}
    

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of UnitTest.cc
//---------------------------------------------------------------------------//
