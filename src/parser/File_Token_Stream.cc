//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file File_Token_Stream.cc
 * \author Kent G. Budge
 * \date Wed Jan 22 15:18:23 MST 2003
 * \brief Definitions of File_Token_Stream methods.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include "File_Token_Stream.hh"

namespace rtt_parser 
{
using std::string;
using std::set;

//-------------------------------------------------------------------------//
/*!
 * \brief Construct an uninitialized File_Token_Stream
 * 
 * Construct a File_Token_Stream that is not yet associated with a file. Use
 * the default Text_Token_Stream user-defined whitespace characters.
 *
 * This function exists primarily to support construction of arrays of
 * File_Token_Streams.  An example of where this might be useful is in serial
 * code that combines output files produced by each processor in a parallel
 * run. 
 */

File_Token_Stream::File_Token_Stream()
{
    Ensure(check_class_invariants());
    Ensure(location() == "<uninitialized>");
}

//-------------------------------------------------------------------------//
/*!
 * \brief Construct a File_Token_Stream from a file.
 * 
 * Construct a File_Token_Stream that derives its text from the
 * specified file. If the file cannot be opened, then \c error()
 * will test true. Use the default Text_Token_Stream user-defined
 * whitespace characters.
 *
 * \param file_name
 * Name of the file from which to extract tokens.
 *
 * \throw std::alloc If there is not enough memory to initialize the queues.
 * \throw rtt_dsxx::assertion If the input stream cannot be opened.
 *
 * \todo Make this constructor failsafe.
 */

File_Token_Stream::File_Token_Stream(string const &file_name)
    :
    filename_(file_name),
    infile_(file_name.c_str())
{
    if( ! infile_ )
    {
	std::ostringstream errmsg;
	errmsg << "Cannot construct File_Token_Stream.\n"
	       << "The file specified could not be found.\n"
	       << "The file requested was: \"" << file_name 
	       << "\"" << std::endl;
	throw rtt_dsxx::assertion( errmsg.str().c_str() );
    }

    Ensure(check_class_invariants());
    Ensure(location() == file_name + ", line 1");
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Construct a File_Token_Stream from a file.
 * 
 * Construct a File_Token_Stream that derives its text from the
 * specified file. If the file cannot be opened, then \c error()
 * will test true. 
 *
 * \param file_name Name of the file from which to extract tokens.
 *
 * \param ws Points to a string containing user-defined whitespace characters.
 */

File_Token_Stream::File_Token_Stream(string const &file_name,
				     set<char> const &ws)
    :
    Text_Token_Stream(ws),
    filename_(file_name),
    infile_(file_name.c_str())
{
    if( ! infile_ )
    {
	std::ostringstream errmsg;
	errmsg << "Cannot construct File_Token_Stream.\n"
	       << "The file specified could not be found.\n"
	       << "The file requested was: \"" << file_name 
	       << "\"" << std::endl;
	throw rtt_dsxx::assertion( errmsg.str().c_str() );
    }

    Ensure(check_class_invariants());
    Ensure(location() == file_name + ", line 1");
    Ensure(Whitespace() == ws);
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Attach the File_Token_Stream to a different file.
 * 
 * \param file_name
 * Name of the file to which to attach the token stream.
 *
 * \throw rtt_dsxx::assertion If the input stream cannot be opened.
 */

void File_Token_Stream::open(string file_name)
{
    infile_.close();
    infile_.clear();
    infile_.open(file_name.c_str());
    filename_ = file_name;
    
    if( ! infile_ )
    {
	std::ostringstream errmsg;
	errmsg << "Cannot open File_Token_Stream.\n"
	       << "The file specified could not be found.\n"
	       << "The file requested was: \"" << file_name 
	       << "\"" << std::endl;
	throw rtt_dsxx::assertion( errmsg.str().c_str() );
    }

    Rewind();

    Ensure(check_class_invariants());
    Ensure(location() == file_name + ", line 1");
}

//-------------------------------------------------------------------------//
/*!
 * \brief Returns a locator string.
 *
 * This function constructs and returns a string of the form
 * "filename, line #" indicating the location at which the last
 * token was parsed.  This is useful for error reporting in parsers.
 *
 * \return A string of the form "filename, line #"
 */

string File_Token_Stream::location() const
{
    Require(check_class_invariants());

    std::ostringstream Result;
    if (filename_ != "")
    {
        Result << filename_ << ", line " << Line();
    }
    else
    {
        Result << "<uninitialized>";
    }
    return Result.str();
}
  
//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Fill the character buffer.
 *
 * This function moves the next character in the file stream into the
 * character buffer.
 */

void File_Token_Stream::fill_character_buffer()
{
    char const c = infile_.get();
    if (infile_.fail())
    {
	character_push_back('\0');
    }
    else
    {
	character_push_back(c);
    }
    
    Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Return error flag.
 *
 * This function may be used to check whether an I/O error has occured,
 * such as failure to open the text file.
 *
 * \return \c true if an error has occured; \c false otherwise.
 */

bool File_Token_Stream::error() const
{
    return infile_.fail();
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Return end of file flag.
 *
 * This function may be used to check whether the end of the text file
 * has been reached.
 *
 * \return \c true if the end of the text file has been reached; \c false
 * otherwise.
 */

bool File_Token_Stream::end() const
{
    return infile_.eof();
}

//-------------------------------------------------------------------------//
/* virtual */
void File_Token_Stream::Report(Token const &token,
                               string const &message)
{
    std::cerr << token.Location() << ": " << message << std::endl;

    Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to the error console stream.
 * This version assumes that the cursor gives the correct error location.
 */

void File_Token_Stream::Report(string const &message)
{
    Token const token = Lookahead();
    std::cerr << token.Location() << ": " << message << std::endl;

    Ensure(check_class_invariants());
}
  
//-------------------------------------------------------------------------//
/*!
 * \brief Rewind the file token stream.
 *
 * This function rewinds the file stream associated with the file token
 * stream and flushes its internal buffers, so that scanning resumes at
 * the beginning of the file stream. The error count is also reset.
 */

void File_Token_Stream::Rewind()
{
    Require(check_class_invariants());

    infile_.clear();    // Must clear the error/end flag bits.
    infile_.seekg(0);

    Text_Token_Stream::Rewind();
    
    Ensure(check_class_invariants());
    Ensure(location() == filename_ + ", line 1");
    Ensure(Error_Count()==0);
}

} // namespace rtt_parser

//---------------------------------------------------------------------------//
//                      end of File_Token_Stream.cc
//---------------------------------------------------------------------------//
