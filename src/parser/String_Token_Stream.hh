//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   parser/String_Token_Stream.hh
 * \author Kent G. Budge
 * \brief  Definition of class String_Token_Stream.
 * \note   Copyright (C) 2006-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef CCS4_String_Token_Stream_HH
#define CCS4_String_Token_Stream_HH

#include <fstream>
#include "Text_Token_Stream.hh"

#if defined(MSVC)
#   pragma warning (push)
#   pragma warning (disable:4251) // warning C4251: 'rtt_parser::Abstract_Class_Parser_Base::keys_' : class 'rtt_parser::Abstract_Class_Parser_Base::c_string_vector' needs to have dll-interface to be used by clients of class 'rtt_parser::Abstract_Class_Parser_Base'
#endif

namespace rtt_parser 
{
using std::string;
using std::set;

//-------------------------------------------------------------------------//
/*! 
 * \brief std::string-based token stream
 *
 * String_Token_Stream is a Text_Token_Stream that obtains its text from a
 * std::string passed to the constructor. The diagnostic output is directed to
 * an internal string that can be retrieved at will.
 */

class DLL_PUBLIC String_Token_Stream : public Text_Token_Stream
{
  public:

    // CREATORS

    //! Construct a String_Token_Stream from a string.
    String_Token_Stream(string const &text);

    //! Construct a String_Token_Stream from a file.
    String_Token_Stream(string const &text,
                        set<char> const &whitespace);

    // MANIPULATORS

    // Return to the start of the string.
    void rewind();

    //! Report a condition.
    virtual void report(Token const & token,
                        string const &message);

    //! Report a condition.
    virtual void report(string const &message);

    // ACCESSORS

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

  private:

    // IMPLEMENTATION

    // DATA

    string text_;       //!< Text to be tokenized
    unsigned pos_;      //!< Cursor position in string

    string messages_;   //!< Collection of diagnostic messages
};

} // rtt_parser

#if defined(MSVC)
#   pragma warning (pop)
#endif

#endif  // CCS4_String_Token_Stream_HH
//---------------------------------------------------------------------------//
// end of String_Token_Stream.hh
//---------------------------------------------------------------------------//
