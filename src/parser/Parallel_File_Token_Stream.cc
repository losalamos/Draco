//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Parallel_File_Token_Stream.cc
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Definitions of Parallel_File_Token_Stream methods.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Parallel_File_Token_Stream.hh"
#include "c4/C4_Functions.hh"
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

namespace rtt_parser {
using namespace std;
using namespace rtt_dsxx;

//-------------------------------------------------------------------------//
/*!
 *
 * Construct an empty Parallel_File_Token_Stream..
 */

Parallel_File_Token_Stream::Parallel_File_Token_Stream()
    : is_io_processor_(rtt_c4::node() == 0),
      // The current implementation always designates processor 0 as the I/O
      // processor.
      at_eof_(true), at_error_(false) {
  Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 *
 * Construct a Parallel_File_Token_Stream that derives its text from the
 * specified file. This Parallel_File_Token_Stream will use the default
 * Text_Token_Stream breaking whitespace characters. If the file cannot be
 * opened, then an exception is thrown.
 *
 * \param file_name Name of the file from which to extract tokens.
 *
 * \throw std::invalid_argument If the file cannot be opened.
 *
 * \todo Make this constructor more failsafe.
 */

Parallel_File_Token_Stream::Parallel_File_Token_Stream(string const &file_name)
    : filename_(file_name), infile_(), is_io_processor_(rtt_c4::node() == 0),
      // The current implementation always designates processor 0 as the I/O
      // processor.
      at_eof_(false), at_error_(false) {
  open_();

  Ensure(check_class_invariants());
  Ensure(location_() == file_name + ", line 1");
}

//-------------------------------------------------------------------------//
/*!
 *
 * Construct a Parallel_File_Token_Stream that derives its text from the
 * specified file. If the file cannot be opened, then \c error()
 * will test true.
 *
 * \param file_name Name of the file from which to extract tokens.
 *
 * \param ws Points to a string containing user-defined whitespace characters.
 *
 * \throw std::bad_alloc If there is not enough memory to initialize the
 * token and character queues.
 * \throw std::invalid_argument If the file cannot be opened.
 *
 * \todo Make this constructor more failsafe.
 */

Parallel_File_Token_Stream::Parallel_File_Token_Stream(string const &file_name,
                                                       set<char> const &ws)
    : Text_Token_Stream(ws), filename_(file_name), infile_(),
      is_io_processor_(rtt_c4::node() == 0), at_eof_(false), at_error_(false) {
  open_();

  Ensure(check_class_invariants());
  Ensure(location_() == file_name + ", line 1");
  Ensure(whitespace() == ws);
}

//-------------------------------------------------------------------------//
/*!
 * Reopen a Parallel_File_Token_Stream with a new file.
 *
 * \param file_name Name of the file from which to extract tokens.
 *
 * \throw std::invalid_argument If the file cannot be opened.
 *
 * \todo Make this constructor more failsafe.
 */

void Parallel_File_Token_Stream::open(string const &file_name) {
  filename_ = file_name;
  infile_.close();
  at_eof_ = false;
  Text_Token_Stream::rewind();

  open_();

  Ensure(check_class_invariants());
  Ensure(location_() == file_name + ", line 1");
}
//---------------------------------------------------------------------------//
/*!
 *
 *
 * \throw std::invalid_argument If the file cannot be opened.
 *
 * \todo Make this function more failsafe.
 */

/* private */
void Parallel_File_Token_Stream::open_() {
  // Create in input stream by opening the specified file on the IO proc.
  at_error_ = false;
  if (is_io_processor_) {
    infile_.open(filename_.c_str());
    if (!infile_) {
      at_error_ = true;
    }
  }
  rtt_c4::broadcast(&at_error_, 1, 0);
  if (at_error_) {
    ostringstream errmsg;
    errmsg << "Cannot construct Parallel_File_Token_Stream.\n"
           << "The file specified could not be found.\n"
           << "The file requested was: \"" << filename_ << "\""
           << " (PE " << rtt_c4::node() << ")\n"
           << "Ensure that the filename includes the full path or "
           << "a relative path from\n"
           << "the binary to the input file (e.g. ../test/deck.inp)."
           << "\n\n"
           << "This error can also occur if you forgot to"
           << "execute the code under\n"
           << "mpirun or prun." << endl;
    throw std::invalid_argument(errmsg.str().c_str());
  }
}

//-------------------------------------------------------------------------//
/*!
 *
 * This function constructs and returns a string of the form "filename, line
 * #" indicating the location at which the last token was parsed.  This is
 * useful for error reporting in parsers.
 *
 * \return A string of the form "filename, line #"
 */

string Parallel_File_Token_Stream::location_() const {
  ostringstream Result;
  Result << filename_ << ", line " << line();
  return Result.str();
}

//-------------------------------------------------------------------------//
/*!
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
 * \return The next character in the text stream.
 *
 * \throw rtt_dsxx::assert If a received message has a length greater than
 * the maximum expected.
 */

void Parallel_File_Token_Stream::fill_character_buffer_() {
  using rtt_c4::broadcast;

  // The first value in the communications buffer will be a status code,
  // which if positive is the number of valid characters ini the
  // buffer. This dictates the maximum size needed for the buffer to be the
  // maximum positive character value, plus one (for the status code
  // itself).
  vector<char> comm_buffer(numeric_limits<signed char>::max() + 1);

  // i points to the current position in the communications buffer. We
  // initialize it to 1 to reserve the first character for the status
  // code.
  unsigned i = 1;
  if (is_io_processor_) {
    // Read up to numeric_limits<signed char>::max() characters from the
    // input file.
    while (i < static_cast<unsigned>(numeric_limits<signed char>::max() + 1)) {
      char const c = infile_.get();
      if (infile_.eof() || infile_.fail())
        break;
      comm_buffer[i++] = c;
    }

    if (i > 1) {
      // If there is an end or error condition, but one or more
      // characters were successfully read prior to encountering the
      // end or error condition, wait to transmit the end or error until
      // the next call to fill_character_buffer.

      // Set the status code to the number of valid characters read.
      comm_buffer[0] = static_cast<char>(i - 1);
    } else if (infile_.eof() && !infile_.bad()) {
      // Normal end of file condition.
      comm_buffer[0] = '\0';
    } else {
      // Something went seriously wrong.
      comm_buffer[0] = static_cast<char>(-1);
    }
  }

  vector<char>::iterator first = comm_buffer.begin();
  vector<char>::iterator last = first + i;

  rtt_c4::broadcast(first, last, first);

  if (comm_buffer[0] == '\0') {
    character_push_back_('\0');
    at_eof_ = true;
  } else if (comm_buffer[0] == static_cast<char>(-1)) {
    character_push_back_('\0');
    at_error_ = true;
  } else {
    // Set i to point to the end of the valid sequence of characters in
    // the communications buffer.
    i = 1 + comm_buffer[0];

    // Make sure this is not past the end of the buffer. This should not
    // be possible unless the data has somewhow become corrupted during
    // transmission.
    if (i > static_cast<unsigned>(numeric_limits<signed char>::max() + 1)) {
      throw runtime_error("interprocessor communications corrupted");
    }

    // Copy the transmitted characters into the local character buffer.

    vector<char>::iterator first = comm_buffer.begin();
    vector<char>::iterator last = first + i;

    for (vector<char>::iterator iter = first + 1; iter != last; ++iter) {
      character_push_back_(*iter);
    }
  }

  Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 *
 * This function may be used to check whether an I/O error has occured, such
 * as failure to open the text file.
 *
 * \return \c true if an error has occured; \c false otherwise.
 */

bool Parallel_File_Token_Stream::error_() const { return at_error_; }

//-------------------------------------------------------------------------//
/*!
 *
 * This function may be used to check whether the end of the text file has
 * been reached.
 *
 * \return \c true if the end of the text file has been reached; \c false
 * otherwise.
 */

bool Parallel_File_Token_Stream::end_() const { return at_eof_; }

//-------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to the error console stream.
 * Only processor 0 writes the message, to avoid many (possibly thousands) of
 * duplicate messages.
 */

void Parallel_File_Token_Stream::report(Token const &token,
                                        string const &message) {
  if (rtt_c4::node() == 0) {
    cerr << token.location() << ": " << message << endl;
  }

  Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to the error console stream.
 * Only processor 0 writes the message, to avoid many (possibly thousands) of
 * duplicate messages.
 *
 * This version assumes that the cursor is the message location.
 */

void Parallel_File_Token_Stream::report(string const &message) {
  Require(check_class_invariants());

  Token token = lookahead();
  // The lookahead must be done on all processors to avoid a potential
  // lockup condition.
  if (rtt_c4::node() == 0) {
    cerr << token.location() << ": " << message << endl;
  }

  Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 *
 * This function rewinds the file stream associated with the file token
 * stream and flushes its internal buffers, so that scanning resumes at
 * the beginning of the file stream.
 */

void Parallel_File_Token_Stream::rewind() {
  if (is_io_processor_) {
    infile_.clear(); // Must clear the error/end flag bits.
    infile_.seekg(0);
  }

  at_eof_ = at_error_ = false;

  Text_Token_Stream::rewind();

  Ensure(check_class_invariants());
  Ensure(location_() == filename_ + ", line 1");
}

//---------------------------------------------------------------------------//
bool Parallel_File_Token_Stream::check_class_invariants() const {
  unsigned iocount = is_io_processor_;
  rtt_c4::global_sum(iocount);
  return iocount == 1;
}

} // namespace rtt_parser

//---------------------------------------------------------------------------//
// end of Parallel_File_Token_Stream.cc
//---------------------------------------------------------------------------//
