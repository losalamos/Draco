//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file String_Token_Stream.hh
 * \author Kent G. Budge
 * \date Wed Jan 22 15:18:23 MST 2003
 * \brief Definition of class String_Token_Stream.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef CCS4_String_Token_Stream_HH
#define CCS4_String_Token_Stream_HH

#include <fstream>
#include "Text_Token_Stream.hh"

namespace rtt_parser 
{
//-------------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief std::string-based token stream
 *
 * String_Token_Stream is a Text_Token_Stream that obtains its text from a
 * std::string passed to the constructor. The diagnostic output is directed to
 * an internal string that can be retrieved at will.
 */

class String_Token_Stream : public Text_Token_Stream
{
  public:

    //! Construct a String_Token_Stream from a string.
    String_Token_Stream(std::string const &text);

    //! Construct a String_Token_Stream from a file.
    String_Token_Stream(std::string const &text,
                        std::set<char> const &whitespace);
    
    void Rewind();

    std::string Get_Messages() const { return messages_; }
    
    virtual void Report(Token const & token,
			      std::string const &message);
    
    virtual void Report(std::string const &message);

    bool check_class_invariants() const;
    
  protected:
    
    virtual std::string location() const;
    
    //! Fill the character buffer.
    virtual void fill_character_buffer();

    virtual bool error() const;
    virtual bool end() const;

  private:

    // IMPLEMENTATION

    // DATA

    std::string text_;  //!< Text to be tokenized
    unsigned pos_;      //!< Cursor position in string

    std::string messages_; //!< Collection of diagnostic messages
};

} // rtt_parser

#endif  // CCS4_String_Token_Stream_HH
//---------------------------------------------------------------------------//
//                      end of String_Token_Stream.hh
//---------------------------------------------------------------------------//
