//----------------------------------*-C++-*----------------------------------//
/*!
 * \file Token.cc
 * \author Kent G. Budge
 * \date 18 Feb 2003
 * \brief Definitions of Token helper functions.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 *
 * This file defines the various token type identification functions
 * used in conjuction with struct Token.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "Token.hh"

namespace rtt_parser 
{

//-----------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Is the argument a valid keyword?
 *
 * \return \c true if the argument points to a string consisting of a
 * sequence of  C++ identifiers separated by single spaces.
 *
 * \pre <code>text!=NULL</code>
 */

bool Is_Keyword_Text(char const *text)
{
    Require(text!=NULL);
    
    char c = *text++;
    while (true){
	if (!isalpha(c) && c!='_') return false;
	while (c=*text++, isalnum(c) || c=='_');
	if (!c) return true;
	if (c!=' ') return false;
	c = *text++;
    }
}

//-----------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Is the argument a valid string constant?
 *
 * \return \c true if the argument points to a string consisting of a
 * single C++ string constant, including the delimiting quotes.
 *
 * \pre <code>text!=NULL</code>
 */

bool Is_String_Text(char const *text)
{
    Require(text!=NULL);
    
    char c = *text++;
    if (c!='"') return false;
    while (true){
	c = *text++;
	if (c==0) return false;
	if (c=='"')
	    return !*text++;
	if (c=='\\'){
	    if (!*text++) return false;
	}
    }
}

//-----------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Is the argument a valid real constant?
 *
 * \return \c true if the argument points to a string consisting of a
 * single C++ floating-point constant.
 *
 * \pre <code>text!=NULL</code>
 */

bool Is_Real_Text(char const *text)
{
    Require(text!=NULL);
    
    char *endtext;
    strtod(text, &endtext);
    return !*endtext;
}

//-----------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Is the argument a valid integer constanta?
 *
 * \return \c true if the argument points to a string consisting of a
 * single C++ integer constant.
 *
 * \pre <code>text!=NULL</code>
 */

bool Is_Integer_Text(char const *text)
{
    Require(text!=NULL);
    
    char *endtext;
    strtol(text, &endtext, 0);
    return !*endtext;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Test equality of two tokens
 * 
 * \param a First token to compare
 * \param b Second token to compare
 *
 * \return \c true if the two tokens are equal.
 */

bool operator==(const Token &a, const Token &b)
{
    return 
	a.Type() == b.Type()  &&  
	a.Text() == b.Text()  &&
	a.Location() == b.Location();
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Check that class invariants are satisfied
 *
 * The invariants all reflect the basic requirement that the token text is
 * consistent with the token type.  For example, if the type is REAL, the
 * text must be a valid C representation of a real number, which can be
 * converted to double using atof.
 * 
 * \return \c true if the invariants are all satisfied; \c false otherwise
 */

bool Token::check_class_invariants() const
{
    return 
	(!(type_==END   || type_==EXIT || type_==ERROR) || text_=="") &&
	(type_!=OTHER   || text_.size()==1)  &&
	(type_!=KEYWORD || Is_Keyword_Text(text_.c_str()))  &&
	(type_!=REAL    || Is_Real_Text(text_.c_str()) )  &&
	(type_!=INTEGER || Is_Integer_Text(text_.c_str()) )  &&
	(type_!=STRING  || Is_String_Text(text_.c_str()) );
}

} // rtt_parser
//---------------------------------------------------------------------------//
//                          end of Token_Stream.cc
//---------------------------------------------------------------------------//
