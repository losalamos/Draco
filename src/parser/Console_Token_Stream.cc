//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file Console_Token_Stream.cc
 * \author Kent G. Budge
 * \date Wed Jan 22 15:18:23 MST 2003
 * \brief Definitions of Console_Token_Stream methods.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
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
#include "Console_Token_Stream.hh"

namespace rtt_parser 
{
using namespace std;
//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Construct a Console_Token_Stream.
 * 
 * Construct a Console_Token_Stream. Use the default Text_Token_Stream
 * user-defined whitespace characters.
 */

Console_Token_Stream::Console_Token_Stream()
{
    Ensure(check_class_invariants());
    Ensure(location() == "input");
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Construct a Console_Token_Stream.
 * 
 * Construct a Console_Token_Stream. 
 *
 * \param ws User-defined whitespace characters.
 */

Console_Token_Stream::Console_Token_Stream(set<char> const &ws)
    : Text_Token_Stream(ws)
{
    Ensure(check_class_invariants());
    Ensure(location() == "input");
    Ensure(Whitespace()==ws);
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Returns a locator string.
 *
 * For a Console_Token_Stream, location is not a terribly  meaningful
 * concept.  So we return "input" as the location, which is true enough.
 *
 * \return The string "input".
 */

std::string Console_Token_Stream::location() const
{
    return "input";
}
  
//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Fill the character buffer.
 *
 * This function moves the next character from cin into the character buffer.
 */

void Console_Token_Stream::fill_character_buffer()
{
    char c = cin.get();
    if (cin.fail())
    {
	character_push_back('\0');
    }
    else
    {
	if (c=='\n')
        {
            c=';';
        }
	character_push_back(c);
    }
    
    Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Return error flag.
 *
 * This function may be used to check whether an I/O error has occured.
 *
 * \return \c true if an error has occured; \c false otherwise.
 */

bool Console_Token_Stream::error() const
{
    return cin.fail();
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Return end of file flag.
 *
 * This function may be used to check whether the user has typed an end of
 * file character (ctrl-D on most Unix systems).
 *
 * \return \c true if an end of file character has been typed; \c false
 * otherwise.
 */

bool Console_Token_Stream::end() const
{
    return cin.eof();
}

//-------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to the error console stream.
 */

void Console_Token_Stream::Report(const Token &token,
                                  const std::string &message)
{
    std::cerr << token.Location() << ": " << message << std::endl;

    Ensure(check_class_invariants());
}

//-------------------------------------------------------------------------//
/*!
 * This function sends a message by writing it to the error console stream.
 * This version assumes that the cursor gives the correct message location.
 */

void Console_Token_Stream::Report(const std::string &message)
{
    Token token = Lookahead();
    std::cerr << token.Location() << ": " << message << std::endl;

    Ensure(check_class_invariants());
}
  
//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 * \brief Rewind the token stream.
 *
 * This function flushes cin and resets the error count.
 */

void Console_Token_Stream::Rewind()
{
    cin.clear();    // Must clear the error/end flag bits.
    cin.seekg(0);

    Text_Token_Stream::Rewind();

    Ensure(check_class_invariants());
    Ensure(cin.rdstate() == 0);
    Ensure(location() == "input");
    Ensure(Error_Count()==0);
}

} // namespace rtt_parser

//---------------------------------------------------------------------------//
//                      end of Console_Token_Stream.cc
//---------------------------------------------------------------------------//
