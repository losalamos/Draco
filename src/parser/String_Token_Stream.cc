//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file String_Token_Stream.cc
 * \author Kent G. Budge
 * \date Wed Jan 22 15:18:23 MST 2003
 * \brief Definitions of String_Token_Stream methods.
 * \note   Copyright @ 2005 The Regents of the University of California.
 *
 * This file defines all methods of class String_Token_Stream.
 *
 * revision history:
 * 0) original
 * 1) kgbudge (03/12/03): Fix indentation. Add additional DBC assertions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>
#include <vector>
#include "c4/global.hh"
#include "c4/ParallelUtils.hh"
#include "String_Token_Stream.hh"

namespace rtt_parser 
{
using namespace std;
using rtt_dsxx::assertion;

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
    : text(text),
      pos(0)
{
    Ensure(check_class_invariants());
    Ensure(Whitespace()==Text_Token_Stream::default_whitespace);
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
 * \param text
 * Text from which to extract tokens.
 * \param ws
 * Points to a string containing user-defined whitespace
 * characters.
 *
 * \post <code>Whitespace()==ws</code>
 */

String_Token_Stream::String_Token_Stream(string const &text,
					 set<char> const &ws)
    : Text_Token_Stream(ws),
      text(text),
      pos(0)
{
    Ensure(check_class_invariants());
    Ensure(Whitespace() == ws);
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Returns a locator string.
 *
 * This function constructs and returns a string of the form
 * "filename, line #" indicating the location at which the last
 * token was parsed.  This is useful for error reporting in parsers.
 *
 * \return A string of the form "filename, line #"
 *
 * \throw bad_alloc If there is not enough memory to store the location string.
 */

string String_Token_Stream::location() const
{
    ostringstream Result;
    Result << "near ";
    for (unsigned i=min(0U,pos-10); i<min(pos+10, unsigned (text.length())); ++i) 
    {
	char const c = text[i];
	Result.put(c);
    }
    return Result.str();
}
  
//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 *
 * \post <code> error() || end() || buffer.size()>old buffer.size()
 */

void String_Token_Stream::fill_character_buffer()
{
    if (pos<text.length())
    {
	character_push_back(text[pos++]);
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
 * String_Token_Stream.
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
 * This function may be used to check whether the end of the text file
 * has been reached.
 *
 * \return \c true if the end of the text string has been reached; \c false
 * otherwise.
 */

bool String_Token_Stream::end() const
{
    return pos>=text.length();
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Report an error to the user.
 *
 * This function reports an error by writing it to an internal string.
 */

void String_Token_Stream::Report_Error(Token const &token,
				       string const &message)
{
    messages += token.Location() + ": " + message + '\n';

    Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Report an error to the user.
 *
 * This function reports an error by writing it to an internal string..
 *
 * This version assumes that the cursor is the error location.
 */

void String_Token_Stream::Report_Error(string const &message)
{
    Token token = Lookahead();
    messages += token.Location() + ": " + message + '\n';

    Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Rewind the file token stream.
 *
 * This function sets pos back to the start of the text string.
 *
 * \post <code>Error_Count() == 0</code>
 */

void String_Token_Stream::Rewind()
{
    pos = 0;

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
    return pos<=text.length();
}

}  // namespace rtt_parser
//---------------------------------------------------------------------------//
//                      end of String_Token_Stream.cc
//---------------------------------------------------------------------------//
