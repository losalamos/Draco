//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file String_Token_Stream.cc
 * \author Kent G. Budge
 * \date Wed Jan 22 15:18:23 MST 2003
 * \brief Definitions of String_Token_Stream methods.
 * \note   Copyright @ 2005 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <vector>
#include "c4/global.hh"
#include "String_Token_Stream.hh"

namespace rtt_parser 
{
using namespace std;
using namespace rtt_dsxx;

//-------------------------------------------------------------------------//
/*!
 * Construct a String_Token_Stream that derives its text from the
 * specified string. Use the default Text_Token_Stream user-defined
 * whitespace characters.
 *
 * \param text
 * Text to be tokenized.
 */

String_Token_Stream::String_Token_Stream(string const &text)
    : text_(text),
      pos_(0)
{
    Ensure(check_class_invariants());
    Ensure(Whitespace()==Text_Token_Stream::default_whitespace);
    Ensure(Get_Messages()=="");
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Construct a String_Token_Stream from a file.
 * 
 * Construct a String_Token_Stream that derives its text from the
 * specified string. 
 *
 * \param text Text from which to extract tokens.
 *
 * \param ws Points to a string containing user-defined whitespace characters.
 */

String_Token_Stream::String_Token_Stream(string const &text,
					 set<char> const &ws)
    : Text_Token_Stream(ws),
      text_(text),
      pos_(0)
{
    Ensure(check_class_invariants());
    Ensure(Whitespace() == ws);
    Ensure(Get_Messages()=="");
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Returns a locator string.
 *
 * This function constructs and returns a string of the form "near <text>"
 * where <text> reproduces the region of text where the last token was parsed.
 * This is useful for error reporting in parsers.
 *
 * \return A string of the form "near <text>"
 *
 * \throw bad_alloc If there is not enough memory to store the location string.
 */

string String_Token_Stream::location() const
{
    ostringstream Result;
    Result << "near ";
    for (unsigned i=min(0U, pos_-10);
         i<min(pos_+10, static_cast<unsigned>(text_.length()));
         ++i) 
    {
	char const c = text_[i];
	Result.put(c);
    }
    return Result.str();
}
  
//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 */

void String_Token_Stream::fill_character_buffer()
{
    if (pos_<text_.length())
    {
	character_push_back(text_[pos_++]);
    }
    else
    {
	character_push_back('\x0');
    }

    Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Return error flag.
 *
 * \return \c false; no error conditions are possible for a
 * String_Token_Stream since it performs no I/O.
 */

bool String_Token_Stream::error() const
{
    return false;
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Return end of file flag.
 *
 * This function may be used to check whether the end of the text string
 * has been reached.
 *
 * \return \c true if the end of the text string has been reached; \c false
 * otherwise.
 */

bool String_Token_Stream::end() const
{
    return pos_>=text_.length();
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Report an error to the user.
 *
 * This function sends a messsage by writing it to an internal string.
 */

void String_Token_Stream::Report(Token const &token,
                                 string const &message)
{
    messages_ += token.Location() + ": " + message + '\n';

    Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Report an error to the user.
 *
 * This function sends a message by writing it to an internal string..
 *
 * This version assumes that the cursor is the error location.
 */

void String_Token_Stream::Report(string const &message)
{
    Token token = Lookahead();
    messages_ += token.Location() + ": " + message + '\n';

    Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Rewind the file token stream.
 *
 * This function sets pos back to the start of the text string.
 */

void String_Token_Stream::Rewind()
{
    pos_ = 0;

    Text_Token_Stream::Rewind();
    
    Ensure(check_class_invariants());
    Ensure(Error_Count()==0);
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Check that the class invariants are satisfied.
 * 
 * \return \c true if the invariants are satisfied; \c false otherwise
 */

bool String_Token_Stream::check_class_invariants() const
{
    return pos_<=text_.length();
}

}  // namespace rtt_parser
//---------------------------------------------------------------------------//
//                      end of String_Token_Stream.cc
//---------------------------------------------------------------------------//
