//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   File_Token_Stream.cc
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Definitions of File_Token_Stream methods.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------------------//

#include "File_Token_Stream.hh"
#include <iostream>
#include <sstream>

namespace rtt_parser {
using namespace std;

//-------------------------------------------------------------------------------------//
/*!
 * Construct a File_Token_Stream that is not yet associated with a file. Use
 * the default Text_Token_Stream user-defined whitespace characters.
 *
 * This function exists primarily to support construction of arrays of
 * File_Token_Streams.  An example of where this might be useful is in serial
 * code that combines output files produced by each processor in a parallel
 * run.
 */

File_Token_Stream::File_Token_Stream(void)
    : filename_(std::string()), infile_() {
  Ensure(check_class_invariants());
  Ensure(location_() == "<uninitialized>");
}

//-------------------------------------------------------------------------------------//
/*!
 * Construct a File_Token_Stream that derives its text from the specified
 * file. If the file cannot be opened, then \c error() will test true. Use the
 * default Text_Token_Stream user-defined whitespace characters.
 *
 * \param file_name
 * Name of the file from which to extract tokens.
 *
 * \throw invalid_argument If the input stream cannot be opened.
 *
 * \todo Make this constructor failsafe.
 */

File_Token_Stream::File_Token_Stream(string const &file_name)
    : filename_(file_name), infile_(file_name.c_str(), std::ios::in) {
  if (!infile_) {
    ostringstream errmsg;
    errmsg << "Cannot construct File_Token_Stream.\n"
           << "The file specified could not be found.\n"
           << "The file requested was: \"" << file_name << "\"" << endl;
    throw invalid_argument(errmsg.str().c_str());
  }

  Ensure(check_class_invariants());
  Ensure(location_() == file_name + ", line 1");
}

//-------------------------------------------------------------------------------------//
/*!
 * Construct a File_Token_Stream that derives its text from the specified file. If the file
 * cannot be opened, then \c error() will test true.
 *
 * \param file_name Name of the file from which to extract tokens.
 *
 * \param ws User-defined whitespace characters.
 *
 * \param no_nonbreaking_ws If true, treat spaces and tabs as breaking whitespace. This has
 * the effect of forcing all keywords to consist of a single identifier.
 *
 * \throw invalid_argument If the input stream cannot be opened.
 */

File_Token_Stream::File_Token_Stream(string const &file_name,
                                     set<char> const &ws,
                                     bool const no_nonbreaking_ws)
    : Text_Token_Stream(ws, no_nonbreaking_ws), filename_(file_name),
      infile_(file_name.c_str()) {
  if (!infile_) {
    ostringstream errmsg;
    errmsg << "Cannot construct File_Token_Stream.\n"
           << "The file specified could not be found.\n"
           << "The file requested was: \"" << file_name << "\"" << endl;
    throw invalid_argument(errmsg.str().c_str());
  }

  Ensure(check_class_invariants());
  Ensure(location_() == file_name + ", line 1");
  Ensure(whitespace() == ws);
  Ensure(this->no_nonbreaking_ws() == no_nonbreaking_ws);
}

//---------------------------------------------------------------------------------------//
/*!
 * \brief Attach the File_Token_Stream to a different file.
 *
 * \param file_name
 * Name of the file to which to attach the token stream.
 *
 * \throw invalid_argument If the input stream cannot be opened.
 */

void File_Token_Stream::open(string const &file_name) {
  infile_.close();
  infile_.clear();
  infile_.open(file_name.c_str());
  filename_ = file_name;

  if (!infile_) {
    ostringstream errmsg;
    errmsg << "Cannot open File_Token_Stream.\n"
           << "The file specified could not be found.\n"
           << "The file requested was: \"" << file_name << "\"" << endl;
    throw invalid_argument(errmsg.str().c_str());
  }

  rewind();

  Ensure(check_class_invariants());
  Ensure(location_() == file_name + ", line 1");
}

//-------------------------------------------------------------------------------------//
/*!
 * This function constructs and returns a string of the form
 * "filename, line #" indicating the location at which the last
 * token was parsed.  This is useful for error reporting in parsers.
 *
 * \return A string of the form "filename, line #"
 */

string File_Token_Stream::location_() const {
  ostringstream Result;
  if (filename_ != "") {
    Result << filename_ << ", line " << line();
  } else {
    Result << "<uninitialized>";
  }
  return Result.str();
}

//-------------------------------------------------------------------------------------//
/*!
 * This function moves the next character in the file stream into the
 * character buffer.
 */

void File_Token_Stream::fill_character_buffer_() {
  char const c = infile_.get();
  if (infile_.fail()) {
    character_push_back_('\0');
  } else {
    character_push_back_(c);
  }

  Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------------------//
/*!
 * This function may be used to check whether an I/O error has occured, such
 * as failure to open the text file.
 *
 * \return \c true if an error has occured; \c false otherwise.
 */

bool File_Token_Stream::error_() const { return infile_.fail(); }

//-------------------------------------------------------------------------------------//
/*!
 * This function may be used to check whether the end of the text file has
 * been reached.
 *
 * \return \c true if the end of the text file has been reached; \c false
 * otherwise.
 */

bool File_Token_Stream::end_() const { return infile_.eof(); }

//-------------------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to the error console stream.
 */
void File_Token_Stream::report(Token const &token, string const &message) {
  cerr << token.location() << ": " << message << endl;

  Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to the error console stream.
 * This version assumes that the cursor gives the correct error location.
 */

void File_Token_Stream::report(string const &message) {
  Token const token = lookahead();
  cerr << token.location() << ": " << message << endl;

  Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------------------//
/*!
 * This function rewinds the file stream associated with the file token
 * stream and flushes its internal buffers, so that scanning resumes at
 * the beginning of the file stream. The error count is also reset.
 */

void File_Token_Stream::rewind() {
  infile_.clear(); // Must clear the error/end flag bits.
  infile_.seekg(0);

  Text_Token_Stream::rewind();

  Ensure(check_class_invariants());
  Ensure(location_() == filename_ + ", line 1");
}

} // namespace rtt_parser

//---------------------------------------------------------------------------------------//
// end of File_Token_Stream.cc
//---------------------------------------------------------------------------------------//
