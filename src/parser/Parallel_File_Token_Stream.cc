//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file Parallel_File_Token_Stream.cc
 * \author Kent G. Budge
 * \date Wed Jan 22 15:18:23 MST 2003
 * \brief Definitions of Parallel_File_Token_Stream methods.
 * \note   Copyright @ 2003 The Regents of the University of California.
 *
 * This file defines all methods of class Parallel_File_Token_Stream.
 *
 * revision history:
 * 0) original
 * 1) kgbudge (03/12/03): Fix indentation. Add additional DBC assertions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <vector>
#include "c4/global.hh"
#include "Parallel_File_Token_Stream.hh"

namespace rtt_parser 
{
using namespace std;
using rtt_dsxx::assertion;

//-------------------------------------------------------------------------//
/*!
 *
 * Construct a Parallel_File_Token_Stream that derives its text from the
 * specified file. This Parallel_File_Token_Stream will use the default
 * Text_Token_Stream breaking whitespace characters. If the file cannot be
 * opened, then an exception is thrown.
 *
 * \param file_name
 * Name of the file from which to extract tokens.
 *
 * \post <code>location() == file_name + ", line 1"</code>
 * \post \c Whitespace()==Text_Token_Stream::default_whitespace
 *
 * \throw assertion If the file cannot be opened.
 *
 * \todo Make this constructor more failsafe.
 */

Parallel_File_Token_Stream::Parallel_File_Token_Stream(string const &file_name)
    : filename(        file_name         ),
      is_io_processor( rtt_c4::node()==0 ),  // The current implementation
					     // always designates processor 0
					     // as the I/O  processor. 
      at_eof(   false ),
      at_error( false )
{
    open();

    Ensure(check_class_invariants());
    Ensure(Whitespace()==Text_Token_Stream::default_whitespace);
    Ensure(location() == file_name + ", line 1");
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Construct a Parallel_File_Token_Stream from a file.
 * 
 * Construct a Parallel_File_Token_Stream that derives its text from the
 * specified file. If the file cannot be opened, then \c error()
 * will test true. 
 *
 * \param file_name
 * Name of the file from which to extract tokens.
 * \param ws
 * Points to a string containing user-defined whitespace
 * characters.
 *
 * \post <code>location() == file_name + ", line 1"</code>
 * \post <code>Whitespace()==ws</code>
 *
 * \throw std::bad_alloc If there is not enough memory to initialize the
 * token and character queues.
 * \throw assertion If the file cannot be opened.
 *
 * \todo Make this constructor more failsafe.
 */

Parallel_File_Token_Stream::Parallel_File_Token_Stream(string const &file_name,
						       set<char> const &ws)
    : Text_Token_Stream(ws),
      filename(file_name),
      is_io_processor(rtt_c4::node()==0),
      at_eof(false),
      at_error(false)
{
    open();

    Ensure(check_class_invariants());
    Ensure(location() == file_name + ", line 1");
    Ensure(Whitespace() == ws);
}

//---------------------------------------------------------------------------//
/*! 
 *
 * \throw assertion If the file cannot be opened.
 *
 * \todo Make this function more failsafe.
 */

void Parallel_File_Token_Stream::open()
{
    // Create in input stream by opening the specified file on the IO proc.
    if( is_io_processor ) 
    {
	infile.open(filename.c_str());
	if (!infile) at_error = true;
    }
    unsigned err_count = at_error;
    rtt_c4::global_sum( err_count );
    if( err_count > 0 )
    {
	ostringstream errmsg;
	errmsg << "Cannot construct Parallel_File_Token_Stream.\n"
	       << "The file specified could not be found.\n"
	       << "The file requested was: \"" << filename 
	       << "\"" << " (PE " << rtt_c4::node() << ")" <<endl;
	throw assertion( errmsg.str().c_str() );
    }
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Returns a locator string.
 *
 * This function constructs and returns a string of the form
 * "filename, line #" indicating the location at which the last
 * token was parsed.  This is useful for error reporting in parsers.
 *
 * \return A string of the form "filename, line #"
 *
 * \throw bad_alloc If there is not enough memory to store the location string.
 */

string Parallel_File_Token_Stream::location() const
{
    ostringstream Result;
    Result << filename << ", line " << Line();
    return Result.str();
}
  
//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 *
 * Only the I/O processor actually reads from the file.  Up to
 * numeric_limits<signed char>::max() characters are read by this processor.
 * The I/O processor then broadcasts a message consisting of a status
 * character and the characters that were read. If the I/O processor has
 * reached the end of the file, the status character is 0. If the I/O
 * processor has encountered some kind of stream error, the status character
 * is set to \c static_cast<char>(-1).  Otherwise the status character is the
 * number of characters to be transmitted.
 *
 * \post <code> error() || end() || buffer.size()>old buffer.size()
 *
 * \return The next character in the text stream.
 *
 * \throw std::bad_alloc If there is not enough memory to expand the
 * character queue.
 * \throw rtt_dsxx::assert If a received message has a length greater than
 * the maximum expected. 
 */

void Parallel_File_Token_Stream::fill_character_buffer()
{
    using rtt_c4::broadcast;
    
    vector<char> comm_buffer(numeric_limits<signed char>::max()+1);
    unsigned i = 1;    // first character is status character
    if (is_io_processor)
    {
	while (i < static_cast<unsigned>(numeric_limits<signed char>::max()+1))
	{
	    char const c = infile.get();
	    if (infile.eof() || infile.fail()) break;
	    comm_buffer[i++] = c;
	}
	
	if (i>1)
	{
	    // If there is an end or error condition, but one or more
	    // characters were successfully read prior to encountering the
	    // end or error condition, wait to transmit the end or error until
	    // the next call to fill_character_buffer.  
	    comm_buffer[0] = static_cast<char>(i-1);
        }   
        else if (infile.eof() && !infile.bad())
	{
	    // Normal end of file condition.
	    comm_buffer[0] = '\0';
	}
	else
	{
	    // Something went seriously wrong.
	    comm_buffer[0] = static_cast<char>(-1);
	}
    }

    vector<char>::iterator first=comm_buffer.begin();
    vector<char>::iterator last=first+i;

    rtt_c4::broadcast(first, last, first);

    if (comm_buffer[0]=='\0')
    {
	character_push_back('\0');
	at_eof = true;
    }
    else if (comm_buffer[0]==static_cast<char>(-1))
    {
	character_push_back('\0');
	at_error = true;
    }
    else
    {
	i = 1 + comm_buffer[0];

	if (i>static_cast<unsigned>(numeric_limits<signed char>::max()+1))
	{
	    throw runtime_error("interprocessor communications corrupted");
	}

        vector<char>::iterator first=comm_buffer.begin();
        vector<char>::iterator last=first+i;

	for (vector<char>::iterator iter = first+1; iter != last; ++iter)
	{
	    character_push_back(*iter);
	}
    }

    rtt_c4::global_barrier();

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

bool Parallel_File_Token_Stream::error() const
{
    return at_error;
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

bool Parallel_File_Token_Stream::end() const
{
    return at_eof;
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Report an error to the user.
 *
 * This function reports an error by writing it to the error console stream.
 * Only processor 0 writes the error, to avoid many (possibly thousands) of
 * duplicate messages.
 */

void Parallel_File_Token_Stream::Report_Error(Token const &token,
					      string const &message)
{
    if (rtt_c4::node()==0)
    {
	cerr << token.Location() << ": " << message << endl;
    }

    Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Report an error to the user.
 *
 * This function reports an error by writing it to the error console stream.
 * Only processor 0 writes the error, to avoid many (possibly thousands) of
 * duplicate messages.
 *
 * This version assumes that the cursor is the error location.
 */

void Parallel_File_Token_Stream::Report_Error(string const &message)
{
    Require(check_class_invariants());

    if (rtt_c4::node()==0)
    {
	Token token = Lookahead();
	cerr << token.Location() << ": " << message << endl;
    }

    Ensure(check_class_invariants());
}
  
//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Rewind the file token stream.
 *
 * This function rewinds the file stream associated with the file token
 * stream and flushes its internal buffers, so that scanning resumes at
 * the beginning of the file stream.
 *
 * \post <code>location() == filename + ", line 1"</code>;
 * \post <code>Error_Count() == 0</code>
 */

void Parallel_File_Token_Stream::Rewind()
{
    Require(check_class_invariants());

    if (is_io_processor)
    {
	infile.clear();    // Must clear the error/end flag bits.
	infile.seekg(0);
    }

    at_eof = at_error = false;

    Text_Token_Stream::Rewind();
    
    Ensure(check_class_invariants());
    Ensure(location() == filename + ", line 1");
    Ensure(Error_Count()==0);
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Check that the class invariants are satisfied.
 * 
 * \return \c true if the invariants are satisfied; \c false otherwise
 */

bool Parallel_File_Token_Stream::check_class_invariants() const
{
    unsigned iocount = is_io_processor;
    rtt_c4::global_sum(iocount);
    return iocount == 1;
}

}  // namespace rtt_parser
//---------------------------------------------------------------------------//
//                      end of Parallel_File_Token_Stream.cc
//---------------------------------------------------------------------------//
