//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   parser/String_Token_Stream.hh
 * \author Kent G. Budge
 * \brief  Definition of class String_Token_Stream.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved.
 */
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

#ifndef rtt_String_Token_Stream_HH
#define rtt_String_Token_Stream_HH

#include "Text_Token_Stream.hh"
#include <fstream>

namespace rtt_parser {
using std::set;
using std::string;

//-------------------------------------------------------------------------//
/*!
 * \brief std::string-based token stream
 *
 * String_Token_Stream is a Text_Token_Stream that obtains its text from a
 * std::string passed to the constructor. The diagnostic output is directed to
 * an internal string that can be retrieved at will.
 */

class DLL_PUBLIC_parser String_Token_Stream : public Text_Token_Stream {
public:
  // CREATORS

  //! Construct a String_Token_Stream from a string.
  String_Token_Stream(string const &text);

  //! Construct a String_Token_Stream from a string.
  String_Token_Stream(string &&text);

  //! Construct a String_Token_Stream from a string.
  String_Token_Stream(string const &text, set<char> const &whitespace,
                      bool no_nonbreaking_ws = false);

  // MANIPULATORS

  // Return to the start of the string.
  void rewind();

  //! Report a condition.
  virtual void report(Token const &token, string const &message);

  //! Report a condition.
  virtual void report(string const &message);

  //! Report a comment.
  virtual void comment(std::string const &message);

  // ACCESSORS

  //! Return the text to be tokenized.
  string const &text() const { return text_; }

  //! Return the accumulated set of messages.
  string messages() const { return messages_; }

  //! Check the class invariant.
  bool check_class_invariants() const;

protected:
  //! Generate a locator string.
  virtual string location_() const;

  virtual void fill_character_buffer_();

  virtual bool error_() const;
  virtual bool end_() const;

  virtual void push_include(std::string &include_file_name);
  virtual void pop_include();

private:
  // IMPLEMENTATION

  // DATA

  string text_;  //!< Text to be tokenized
  unsigned pos_; //!< Cursor position in string

  string messages_; //!< Collection of diagnostic messages
};

} // namespace rtt_parser

#endif // rtt_String_Token_Stream_HH

//----------------------------------------------------------------------------//
// end of String_Token_Stream.hh
//----------------------------------------------------------------------------//
