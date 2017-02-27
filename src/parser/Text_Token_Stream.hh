//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/Text_Token_Stream.hh
 * \author Kent G. Budge
 * \brief  Definition of the Text_Token_Stream class.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef CCS4_Text_Token_Stream_HH
#define CCS4_Text_Token_Stream_HH

#include "Token_Stream.hh"
#include <set>

namespace rtt_parser {
//-------------------------------------------------------------------------//
/*!
 * \brief Abstract text-based token stream for simple parsers.
 *
 * A Text_Token_Stream constructs its stream of tokens by scanning a stream of
 * text characters, supplied through the \c getc method of its child classes.
 *
 * C and C++ comments are treated as breaking whitespace.
 *
 * Null characters are not permitted in the character stream.  They are used
 * internally to indicate the end of file or an error condition.
 */

class DLL_PUBLIC_parser Text_Token_Stream : public Token_Stream {
public:
  // ACCESSORS

  //! Return the current line in the text stream.
  unsigned line() const {
    Ensure(line_ > 0);
    return line_;
  }

  bool no_nonbreaking_ws() const { return no_nonbreaking_ws_; }

  //! Return the current set of whitespace characters.
  std::set<char> const &whitespace() const { return whitespace_; }

  //! Check the class invariants.
  bool check_class_invariants() const;

  // MANIPULATORS

  virtual void rewind() = 0;

  // SERVICES

  //! Does the Token_Stream consider \c c to be whitespace?
  bool is_whitespace(char c) const;

  //! Does the Token_Stream consider <i>c</i> to be nonbreaking
  bool is_nb_whitespace(char c) const;

  // CONST DATA

  //! The default whitespace definition
  static std::set<char> const default_whitespace;

protected:
  // IMPLEMENTATION

  //! Construct a Text_Token_Stream.
  Text_Token_Stream();

  //! Construct a Text_Token_Stream.
  Text_Token_Stream(std::set<char> const &, bool no_nonbreaking_ws = false);

  //! Scan the next token.
  virtual Token fill_();

  //! Push a character onto the back of the character queue.
  void character_push_back_(char c);

  //! Move one or more characters from the text stream into the character
  //! buffer.
  virtual void fill_character_buffer_() = 0;

  virtual bool error_() const = 0;
  virtual bool end_(void) const = 0;
  virtual std::string location_(void) const = 0;

  //! Pop a character off the internal buffer.
  char pop_char_(void);
  //! Peek ahead at the internal buffer.
  char peek_(unsigned pos = 0);

  //! Skip any whitespace at the cursor position.
  void eat_whitespace_(void);

  // The following scan_ functions are for numeric scanning.  The names
  // reflect the context-free grammar given by Stroustrup in appendix A
  // of _The C++ Programming Language_.  However, we do not presently
  // recognize type suffixes on either integers or floats.
  unsigned scan_floating_literal_(void);
  unsigned scan_digit_sequence_(unsigned &);
  unsigned scan_exponent_part_(unsigned &);
  unsigned scan_fractional_constant_(unsigned &);

  unsigned scan_integer_literal_();
  unsigned scan_decimal_literal_(unsigned &);
  unsigned scan_hexadecimal_literal_(unsigned &);
  unsigned scan_octal_literal_(unsigned &);

private:
  // IMPLEMENTATION

  // DATA

  std::deque<char> buffer_;
  //!< Character buffer. Refilled as needed using fill_character_buffer_()

  std::set<char> whitespace_;
  //!< The whitespace character list

  unsigned line_;
  //!< Current line in input file.

  bool no_nonbreaking_ws_;
  //!< Treat all whitespace as breaking whitespace.
};

} // namespace rtt_parser

#endif // CCS4_Text_Token_Stream_HH
//--------------------------------------------------------------------//
// end of Text_Token_Stream.hh
//--------------------------------------------------------------------//
