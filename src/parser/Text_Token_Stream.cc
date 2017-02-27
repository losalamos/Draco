//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   Text_Token_Stream.cc
 * \author Kent G. Budge
 * \brief  Contains definitions of all Text_Token_Stream member functions.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//----------------------------------------------------------------------------//
// $Id$
//----------------------------------------------------------------------------//

#include "Text_Token_Stream.hh"
#include <cstring>
#include <ctype.h>
#include <string>

namespace rtt_parser {
using namespace std;

//----------------------------------------------------------------------------//
char const default_ws_string[] = "=:;,";

set<char> const Text_Token_Stream::default_whitespace(
    default_ws_string, default_ws_string + sizeof(default_ws_string));

//----------------------------------------------------------------------------//
/*!
 * \brief Constructs a Text_Token_Stream with the specified set of breaking
 * whitespace characters.
 *
 * \param ws String containing the user-defined whitespace characters for this
 * Text_Token_Stream.
 *
 * \param no_nonbreaking_ws If true, treat spaces and tabs as breaking
 * whitespace. This has the effect of forcing all keywords to consist of a
 * single identifier.
 *
 * Whitespace characters are classified as breaking or nonbreaking whitespace.
 * Nonbreaking whitespace separates non-keyword tokens and identifiers within a
 * keyword but has no other significance. Breaking whitespace is similar to
 * nonbreaking whitespace except that it always separates tokens; thus, two
 * identifiers separated by breaking whitespace are considered to belong to
 * separate keywords.
 *
 * Nonbreaking whitespace characters are the space and horizontal tab
 * characters.
 *
 * Breaking whitespace characters include all other characters for which the
 * standard C library function <CODE>isspace(char)</CODE> returns a nonzero
 * value, plus additional characters defined as nonbreaking whitespace by the
 * client of the Token_Stream. In particular, a newline character is always
 * breaking whitespace.
 *
 * Whitespace is stripped from the beginning and end of every token, and the
 * nonbreaking whitespace separating each identifier within a keyword is
 * replaced by a single space character.
 */

Text_Token_Stream::Text_Token_Stream(set<char> const &ws,
                                     bool const no_nonbreaking_ws)
    : buffer_(), whitespace_(ws), line_(1),
      no_nonbreaking_ws_(no_nonbreaking_ws) {
  Ensure(check_class_invariants());
  Ensure(ws == whitespace());
  Ensure(line() == 1);
  Ensure(this->no_nonbreaking_ws() == no_nonbreaking_ws);
}

//----------------------------------------------------------------------------//
/*!
 * \brief Constructs a Text_Token_Stream with the default set of breaking
 * whitespace characters. See the previous constructor documentation for a
 * discussion of how whitespace is defined and handled.
 *
 * The default whitespace characters are contained in the set
 * \c Text_Token_Stream::default_whitespace.
 */

Text_Token_Stream::Text_Token_Stream(void)
    : buffer_(), whitespace_(default_whitespace), line_(1),
      no_nonbreaking_ws_(false) {
  Ensure(check_class_invariants());
  Ensure(whitespace() == default_whitespace);
  Ensure(line() == 1);
}

//----------------------------------------------------------------------------//
/*!
 * \brief Scan the next token from the character stream. The character stream is
 * accessed via the fill_character_buffer, error, and end functions, which are
 * pure virtual functions.
 */

Token Text_Token_Stream::fill_() {
  eat_whitespace_();

  char c = peek_(); // Character at the current cursor position

  string token_location = location_();

  Token returnValue(END, token_location);

  if (c == '\0') {
    // Sentinel value for error or end of file.
    if (end_()) {
      Ensure(check_class_invariants());
      returnValue = Token(EXIT, token_location);
    } else {
      Ensure(check_class_invariants());
      returnValue = Token(rtt_parser::ERROR, token_location);
    }
  } else {
    if (isalpha(c) || c == '_')
    // Beginning of a keyword or END token
    {
      string text(1, c);
      pop_char_();
      c = peek_();
      do {
        // Scan a C identifier.
        while (isalnum(c) || c == '_') {
          text += c;
          pop_char_();
          c = peek_();
        }
        if (!no_nonbreaking_ws_) {
          // Replace any nonbreaking whitespace after the identifier
          // with a single space, but ONLY if the identifier is
          // followed by another identifer.
          while (is_nb_whitespace(c)) {
            pop_char_();
            c = peek_();
          }
          if (isalpha(c) || c == '_')
            text += ' ';
        }
      } while (isalpha(c) || c == '_');

      if (text == "end") {
        Ensure(check_class_invariants());
        return Token(END, token_location);
      } else {
        Ensure(check_class_invariants());
        return Token(KEYWORD, text, token_location);
      }
    } else if (isdigit(c) || c == '.') {
      // A number of some kind.  Note that an initial sign ('+' or '-')
      // is tokenized independently, because it could be interpreted as
      // a binary operator in arithmetic expressions.  It is up to the
      // parser to decide if this is the correct interpretation.
      string text;
      unsigned const float_length = scan_floating_literal_();
      unsigned const int_length = scan_integer_literal_();
      if (float_length > int_length) {
        for (unsigned i = 0; i < float_length; i++) {
          c = pop_char_();
          text += c;
        }
        Ensure(check_class_invariants());
        return Token(REAL, text, token_location);
      } else if (int_length > 0) {
        for (unsigned i = 0; i < int_length; i++) {
          char c = pop_char_();
          text += c;
        }
        Ensure(check_class_invariants());
        return Token(INTEGER, text, token_location);
      } else {
        Check(c == '.');
        pop_char_();
        Ensure(check_class_invariants());
        return Token('.', token_location);
      }
    } else if (c == '"')
    // Manifest string
    {
      string text(1, c);
      pop_char_();
      c = peek_();
      for (;;) {
        while (c != '"' && c != '\\' && c != '\n' && !end_() && !error_()) {
          text += c;
          pop_char_();
          c = peek_();
        }
        if (c == '"')
          break;
        if (c == '\\') {
          text += c;
          pop_char_();
          c = pop_char_();
          text += c;
          c = peek_();
        } else {
          if (end_() || error_()) {
            report_syntax_error(Token(EXIT, token_location),
                                "unexpected end of file; "
                                "did you forget a closing quote?");
          } else {
            Check(c == '\n');
            report_syntax_error(Token(EXIT, token_location),
                                "unexpected end of line; "
                                "did you forget a closing quote?");
          }
        }
      }
      text += '"';
      pop_char_();
      Ensure(check_class_invariants());
      return Token(STRING, text, token_location);
    } else if (c == '<')
    // Multicharacter OTHER
    {
      pop_char_();
      if (peek_() == '=') {
        pop_char_();
        Ensure(check_class_invariants());
        return Token(OTHER, "<=", token_location);
      } else {
        Ensure(check_class_invariants());
        return Token(c, token_location);
      }
    } else if (c == '>')
    // Multicharacter OTHER
    {
      pop_char_();
      if (peek_() == '=') {
        pop_char_();
        Ensure(check_class_invariants());
        return Token(OTHER, ">=", token_location);
      } else {
        Ensure(check_class_invariants());
        return Token(c, token_location);
      }
    } else if (c == '&')
    // Multicharacter OTHER
    {
      pop_char_();
      if (peek_() == '&') {
        pop_char_();
        Ensure(check_class_invariants());
        return Token(OTHER, "&&", token_location);
      } else {
        Ensure(check_class_invariants());
        return Token(c, token_location);
      }
    } else if (c == '|')
    // Multicharacter OTHER
    {
      pop_char_();
      if (peek_() == '|') {
        pop_char_();
        Ensure(check_class_invariants());
        return Token(OTHER, "||", token_location);
      } else {
        Ensure(check_class_invariants());
        return Token(c, token_location);
      }
    } else {
      // OTHER
      pop_char_();
      Ensure(check_class_invariants());
      return Token(c, token_location);
    }
  }
  return returnValue;
}

//----------------------------------------------------------------------------//
/*!
 * \brief This function searches for the argument character in its internal list
 * of whitespace characters.
 *
 * \param c
 * Character to be checked against the whitespace list.
 *
 * \return \c true if and only if the character is found in the internal
 * whitespace list.
 */

bool Text_Token_Stream::is_whitespace(char const c) const {
  return isspace(c) || whitespace_.count(c);
}

//----------------------------------------------------------------------------//
/*!
 * \brief This function searches for the argument character in the
 * Token_Stream's internal list of nonbreaking whitespace characters.
 *
 * \param c Character to be checked against the nonbreaking whitespace list.
 *
 * \return \c true if and only if the character is found in the internal
 * nonbreaking whitespace list, and is \e not found in the breaking whitespace
 * list..
 */

bool Text_Token_Stream::is_nb_whitespace(char const c) const {
  return !whitespace_.count(c) && (c == ' ' || c == '\t');
}

//----------------------------------------------------------------------------//
/*!
 * An internal buffer is used to implement unlimited lookahead, necessary for
 * scanning numbers (which have a quite complex regular expression.)  This
 * function pops a character off the top of the internal buffer, using
 * fill_character_buffer() if necessary to ensure that there is at least one
 * character in the buffer.  If the next character is a carriage return, the
 * line count is incremented.
 *
 * \return The next character in the buffer.
 */

char Text_Token_Stream::pop_char_() {
  Remember(unsigned const old_line = line_);

  char const Result = peek_();
  buffer_.pop_front();
  if (Result == '\n')
    line_++;

  Ensure(check_class_invariants());
  Ensure((Result == '\n' && line_ == old_line + 1) || line_ == old_line);
  return Result;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Try to scan a floating literal.
 *
 * \return length of scanned literal; 0 if no literal could be scanned.
 */

unsigned Text_Token_Stream::scan_floating_literal_() {
  unsigned pos = 0;
  if (scan_fractional_constant_(pos) > 0) {
    scan_exponent_part_(pos);
    return pos;
  } else if (scan_digit_sequence_(pos)) {
    if (scan_exponent_part_(pos) == 0)
      return 0;
    return pos;
  } else {
    return 0;
  }
}

//----------------------------------------------------------------------------//
/*!
 * \brief Try to scan a digit sequence.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_digit_sequence_(unsigned &pos) {
  unsigned const old_pos = pos;
  while (isdigit(peek_(pos)))
    pos++;
  return pos - old_pos;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Try to scan an exponent part.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_exponent_part_(unsigned &pos) {
  unsigned const old_pos = pos;
  char c = peek_(pos);
  if (c == 'e' || c == 'E') {
    pos++;
    c = peek_(pos);
    if (c == '-' || c == '+')
      pos++;
    if (!scan_digit_sequence_(pos)) {
      pos = old_pos;
    }
  }
  return pos - old_pos;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Try to scan a fractional constant.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_fractional_constant_(unsigned &pos) {
  unsigned const old_pos = pos;
  if (scan_digit_sequence_(pos) > 0) {
    if (peek_(pos) != '.') {
      pos = old_pos;
    } else {
      pos++;
      scan_digit_sequence_(pos);
    }
  } else if (peek_(pos) == '.') {
    pos++;
    if (scan_digit_sequence_(pos) == 0) {
      pos = old_pos;
    }
  }
  return pos - old_pos;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Try to scan an integer literal.
 *
 * \return length of scanned literal; 0 if no literal could be scanned.
 */

unsigned Text_Token_Stream::scan_integer_literal_() {
  unsigned pos = 0;
  if (scan_decimal_literal_(pos) > 0) {
  } else if (scan_hexadecimal_literal_(pos)) {
  } else if (scan_octal_literal_(pos)) {
  } else {
    return 0;
  }
  return pos;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Try to scan decimal literal.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_decimal_literal_(unsigned &pos) {
  unsigned const old_pos = pos;
  char c = peek_(pos);
  if (isdigit(c) && c != '0') {
    while (isdigit(c)) {
      pos++;
      c = peek_(pos);
    }
  }
  return pos - old_pos;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Try to scan hexadecimal literal.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_hexadecimal_literal_(unsigned &pos) {
  unsigned old_pos = pos;
  if (peek_(pos) == '0') {
    pos++;
    char c = peek_(pos);
    if (c == 'x' || c == 'X') {
      pos++;
      c = peek_(pos);
      if (!isxdigit(c)) {
        pos = old_pos;
      } else {
        while (isxdigit(c)) {
          pos++;
          c = peek_(pos);
        }
      }
    } else {
      pos = old_pos;
    }
  }
  return pos - old_pos;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Try to scan octal literal.
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_octal_literal_(unsigned &pos) {
  unsigned old_pos = pos;
  char c = peek_(pos);
  if (c == '0') {
    while (isdigit(c) && c != '8' && c != '9') {
      pos++;
      c = peek_(pos);
    }
  }
  return pos - old_pos;
}

//----------------------------------------------------------------------------//
/*!
 * An internal buffer is used to implement unlimited lookahead, necessary for
 * scanning numbers (which have a quite complex regular expression.)  This
 * function peeks ahead the specified number of places in the buffer, using
 * fill_character_buffer() if necessary to ensure that there is a sufficient
 * number of characters in the buffer.
 *
 * \param pos
 * Position at which to peek.
 *
 * \return The character at \c buffer[pos].
 */

char Text_Token_Stream::peek_(unsigned const pos) {
  while (buffer_.size() <= pos) {
    fill_character_buffer_();
  }

  Ensure(check_class_invariants());
  return buffer_[pos];
}

//----------------------------------------------------------------------------//
/*!
 * This function flushes the Text_Token_Stream's internal buffers, so that
 * scanning resumes at the beginning of the file stream.  It is normally called
 * by its overriding version in children of Text_Token_Stream.
 */

void Text_Token_Stream::rewind() {
  buffer_.clear();
  line_ = 1;

  Token_Stream::rewind();

  Ensure(check_class_invariants());
  Ensure(error_count() == 0);
}

//----------------------------------------------------------------------------//
bool Text_Token_Stream::check_class_invariants() const { return line_ > 0; }

//----------------------------------------------------------------------------//
/*!
 * This function skips past any whitespace present at the cursor position,
 * leaving the cursor at the first non-whitespace character following the
 * initial cursor position.
 */
/* private */
void Text_Token_Stream::eat_whitespace_() {
  for (;;) {
    // Scan whitespace
    char c = peek_();
    while (is_whitespace(c) && c != '\0') {
      pop_char_();
      c = peek_();
    }

    // Check for a comment
    if (c == '/') {
      if (peek_(1) == '/') {
        // C++ comment
        while (c != '\n' && !error_() && !end_()) {
          pop_char_();
          c = peek_();
        }
      } else if (peek_(1) == '*') {
        pop_char_(); // pop the '/'
        pop_char_(); // pop the '*'
        while ((peek_(0) != '*' || peek_(1) != '/') && !error_() && !end_()) {
          pop_char_();
        }
        pop_char_(); // pop the '*'
        pop_char_(); // pop the '/'
      } else {
        break;
      }
    } else {
      break;
    }
  }
  // private member function -- no invariant check
}

//----------------------------------------------------------------------------//
/*!
 * \func Text_Token_Stream::location
 *
 * This function returns a location string whose exact format is
 * stream-specific.  For example, for a token stream that scans tokens from a
 * text file, this could be a string of the form "filename, line #".
 *
 * \return A string describing the location from which the Text_Token_Stream is
 * currently scanning tokens.
 */

//----------------------------------------------------------------------------//
/*!
 * \param c Character to be pushed onto the back of the character queue.
 */

void Text_Token_Stream::character_push_back_(char const c) {
  buffer_.push_back(c);

  Ensure(check_class_invariants());
}

} // end namespace rtt_parser

//----------------------------------------------------------------------------//
// end of Text_Token_Stream.cc
//----------------------------------------------------------------------------//
