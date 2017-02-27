//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   String_Token_Stream.cc
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  Definitions of String_Token_Stream methods.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------------------//

#include "String_Token_Stream.hh"
#include "c4/C4_Functions.hh"
#include <sstream>

namespace rtt_parser {
using namespace std;
using namespace rtt_dsxx;

//-------------------------------------------------------------------------------------//
/*!
 * Construct a String_Token_Stream that derives its text from the specified
 * string. Use the default Text_Token_Stream user-defined whitespace
 * characters.
 *
 * \param text
 * Text to be tokenized.
 */

String_Token_Stream::String_Token_Stream(string const &text)
    : text_(text), pos_(0), messages_(std::string()) {
  Ensure(check_class_invariants());
  Ensure(whitespace() == Text_Token_Stream::default_whitespace);
  Ensure(this->text() == text);
  Ensure(messages() == "");
}

//-------------------------------------------------------------------------------------//
/*!
 * Construct a String_Token_Stream that derives its text from the specified
 * string.
 *
 * \param text Text from which to extract tokens.
 *
 * \param ws Points to a string containing user-defined whitespace characters.
 *
 * \param no_nonbreaking_whitespace Causes spaces and tabs to be treated as breaking
 * whitespace. This has the effect of forcing all keywords to consist of a single
 * identifier.
 */

String_Token_Stream::String_Token_Stream(string const &text,
                                         set<char> const &ws,
                                         bool const no_nonbreaking_ws)
    : Text_Token_Stream(ws, no_nonbreaking_ws), text_(text), pos_(0),
      messages_(std::string()) {
  Ensure(check_class_invariants());
  Ensure(whitespace() == ws);
  Ensure(this->text() == text);
  Ensure(messages() == "");
}

//-------------------------------------------------------------------------------------//
/*!
 * This function constructs and returns a string of the form "near <text>"
 * where <text> reproduces the region of text where the last token was parsed.
 * This is useful for error reporting in parsers.
 *
 * \return A string of the form "near <text>"
 */

string String_Token_Stream::location_() const {
  ostringstream Result;
  Result << "near\n";

  // search backwards four endlines
  unsigned begin;
  unsigned count = 0;
  for (begin = pos_; begin > 0; --begin) {
    if (text_[begin] == '\n') {
      if (++count == 4)
        break;
    }
  }
  unsigned const end = text_.size();
  for (unsigned i = begin; i < end; ++i) {
    char const c = text_[i];
    if (i >= pos_ && c == '\n') {
      break;
    }
    Result.put(c);
  }
  Result.put('\n');
  return Result.str();
}

//-------------------------------------------------------------------------------------//
void String_Token_Stream::fill_character_buffer_() {
  if (pos_ < text_.length()) {
    character_push_back_(text_[pos_++]);
  } else {
    character_push_back_('\x0');
  }

  Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------------------//
bool String_Token_Stream::error_() const { return false; }

//-------------------------------------------------------------------------------------//
bool String_Token_Stream::end_() const { return pos_ >= text_.length(); }

//-------------------------------------------------------------------------------------//
/*!
 * This function sends a messsage by writing it to an internal string.
 */

void String_Token_Stream::report(Token const &token, string const &message) {
  messages_ += token.location() + "\n" + message + '\n';

  Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to an internal string..
 *
 * This version assumes that the cursor is the error location.
 */

void String_Token_Stream::report(string const &message) {
  Token token = lookahead();
  messages_ += token.location() + "\n" + message + '\n';

  Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------------------//
void String_Token_Stream::rewind() {
  pos_ = 0;

  Text_Token_Stream::rewind();

  Ensure(check_class_invariants());
}

//---------------------------------------------------------------------------------------//
bool String_Token_Stream::check_class_invariants() const {
  return pos_ <= text_.length();
}

} // namespace rtt_parser

//---------------------------------------------------------------------------------------//
// end of String_Token_Stream.cc
//---------------------------------------------------------------------------------------//
