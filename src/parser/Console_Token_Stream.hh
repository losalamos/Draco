//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   parser/Console_Token_Stream.hh
 * \author Kent G. Budge
 * \brief  Definition of class Console_Token_Stream.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved.
 */
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

#ifndef rtt_Console_Token_Stream_HH
#define rtt_Console_Token_Stream_HH

#include "Text_Token_Stream.hh"
#include <fstream>

namespace rtt_parser {
using std::set;
using std::string;

//-------------------------------------------------------------------------//
/*!
 * \brief Console-based token stream
 *
 * Console_Token_Stream represents a text token stream that derives its text
 * stream from the standard console input stream \c cin.  It reports errors
 * to the standard  console error stream \c cerr.
 *
 * This stream also differs from other streams in that the endline character
 * is converted to the semicolon character.  Parsers for use with console
 * streams typically treat the semicolon as an "end of statement" character
 * by specifying that it is NOT a whitespace character and looking for it as
 * a statement terminator.
 *
 * \note This class is an experimental concept and <i> should not be used in
 * production codes </i>.  In particular, the class cannot readily be tested
 * under our current unit testing system, since it is inherently interactive.
 */

class DLL_PUBLIC_parser Console_Token_Stream : public Text_Token_Stream {
public:
  // CREATORS

  //! Construct a Console_Token_Stream.
  Console_Token_Stream();

  //! Construct a Console_Token_Stream.
  explicit Console_Token_Stream(set<char> const &whitespace,
                                bool no_nonbreaking_ws = false);

  // MANIPULATORS

  void rewind();

  virtual void report(const Token &token, const string &message);

  virtual void report(const string &message);

  virtual void comment(std::string const &message);

protected:
  // IMPLEMENTATION

  //! Return a locator string.
  virtual string location_() const;

  virtual void fill_character_buffer_();
  virtual bool error_() const;
  virtual bool end_() const;

  virtual void push_include(std::string &include_file_name);
  virtual void pop_include();
};

} // namespace rtt_parser

#endif // rtt_Console_Token_Stream_HH
//----------------------------------------------------------------------------//
// end of Console_Token_Stream.hh
//----------------------------------------------------------------------------//
