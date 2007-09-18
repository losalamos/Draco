//----------------------------------*-C++-*----------------------------------//
/*!
 * \file Text_Token_Stream.cc
 * \author Kent G. Budge
 * \date Wed Jan 22 14:57:49 MST 2003
 * \brief Contains definitions of all Text_Token_Stream member functions.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <ctype.h>
#include <string>
#include <string.h>
#include "Text_Token_Stream.hh"

namespace rtt_parser 
{
using std::string;
using std::set;

//---------------------------------------------------------------------------//
char const default_ws_string[] = "=:;,";

set<char> const
Text_Token_Stream::
default_whitespace(default_ws_string, 
		   default_ws_string+sizeof(default_ws_string));


//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 *
 * Constructs a Text_Token_Stream with the specified set of breaking
 * whitespace characters.
 *
 * \param ws String containing the user-defined whitespace characters for this
 * Text_Token_Stream.
 *
 * Whitespace characters are classified as breaking or nonbreaking whitespace.
 * Nonbreaking whitespace separates non-keyword tokens and identifiers within
 * a keyword but has no other significance. Breaking whitespace is similar to
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
 * client of the Token_Stream.
 * 
 * Whitespace is stripped from the beginning and end of every token, and the
 * nonbreaking whitespace separating each identifier within a keyword is
 * replaced by a single space character.
 *
 * \post <CODE>Whitespace()==ws</code>
 * \post \c Line()==1
 */

Text_Token_Stream::Text_Token_Stream(set<char> const &ws)
    : whitespace_(ws), line_(1)
{
    Ensure(check_class_invariants());
    Ensure(ws == Whitespace());
    Ensure(Line()==1);
}

//-------------------------------------------------------------------------//
/*!
 * \author Kent G. Budge
 * \date Wed Jan 22 15:35:42 MST 2003
 *
 * Constructs a Text_Token_Stream with the default set of breaking whitespace
 * characters.
 *
 * Whitespace characters are classified as breaking or nonbreaking whitespace.
 * Nonbreaking whitespace separates non-keyword tokens and identifiers within
 * a keyword but has no other significance. Breaking whitespace is similar to
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
 * client of the Token_Stream.
 * 
 * Whitespace is stripped from the beginning and end of every token, and the
 * nonbreaking whitespace separating each identifier within a keyword is
 * replaced by a single space character.
 *
 * By default, the equals sign '=', the colon and semicolon ':' and ';', and
 * the comma ',' are defined as breaking whitespace.  Thus, by default, the
 * sequence <CODE>MY</code> <code>COLOR = GREEN</CODE> is parsed as the
 * keyword token <CODE>MY</code> <code>COLOR</CODE> followed by the keyword
 * token <CODE>GREEN</CODE>.
 *
 * \post \c Line()==1
 * \post \c Whitespace()==default_whitespace
 */

Text_Token_Stream::Text_Token_Stream()
    : whitespace_(default_whitespace), line_(1)
{
    Ensure(check_class_invariants());
    Ensure(Whitespace()==default_whitespace);
    Ensure(Line()==1);
}

//-------------------------------------------------------------------------//
/*!
 * The next token is scanned from the character stream, which is accessed via
 * the fill_character_buffer, error, and end functions.
 *
 * \post <code>size()==old size()+1</code>
 */

Token Text_Token_Stream::fill()
{
    eat_whitespace();

    char c = peek();    // Character at the current cursor position

    string token_location = location();
    
    if (c=='\0')
        // Sentinel value for error or end of file.
    {
	if (end())
	{
	    Ensure(check_class_invariants());
	    return Token(EXIT, token_location);
	}
	else
	{
	    Ensure(check_class_invariants());
	    return Token(ERROR, token_location);
	}
    }
    else 
    {
	if (isalpha(c) || c=='_')
	    // Beginning of a keyword or END token
	{
	    string text(1, c);
	    pop_char();
	    c = peek();
	    do 
	    {
		// Scan a C identifier.
		while (isalnum(c) || c=='_')
		{
		    text += c;
		    pop_char();
		    c = peek();
		}
		// Replace any nonbreaking whitespace after the identifier
		// with a single space, but ONLY if the identifier is
		// followed by another identifer.
		while (is_nb_whitespace(c))
		{
		    pop_char();
		    c = peek();
		}
		if (isalpha(c) || c=='_') text += ' ';
	    } 
	    while (isalpha(c) || c=='_');
	    
	    if (text=="end")
	    {
		Ensure(check_class_invariants());
		return Token(END, token_location);
	    }
	    else
	    {
		Ensure(check_class_invariants());
		return Token(KEYWORD, text, token_location);
	    }
	}
	else if (isdigit(c) || c=='.') 
	    // A number of some kind.  Note that an initial sign ('+' or '-')
	    // is tokenized independently, because it could be interpreted as
	    // a binary operator in arithmetic expressions.  It is up to the
	    // parser to decide if this is the correct interpretation.
	{
	    string text;
	    unsigned const float_length = scan_floating_literal();
	    unsigned const int_length = scan_integer_literal();
	    if (float_length>int_length)
	    {
		for (unsigned i=0; i<float_length; i++)
		{
		    c = pop_char();
		    text += c;
		}
		Ensure(check_class_invariants());
		return Token(REAL, text, token_location);
	    }
	    else if (int_length>0)
	    {
		for (unsigned i=0; i<int_length; i++)
		{
		    char c = pop_char();
		    text += c;
		}
		Ensure(check_class_invariants());
		return Token(INTEGER, text, token_location);
	    }
	    else 
	    {
                Check(c=='.');
                pop_char();
                Ensure(check_class_invariants());
                return Token('.', token_location);
	    }
	}
	else if (c=='"')
	    // Manifest string
	{
	    string text(1, c);
	    pop_char();
	    c = peek();
	    for (;;)
	    {
		while (c!='"' && c!='\\' && c!='\n' && !end() && !error())
		{
		    text += c;
		    pop_char();
		    c = peek();
		}
		if (c=='"') break;
		if (c=='\\')
		{
		    text += c;
		    pop_char();
		    c = pop_char();
		    text += c;
		    c = peek();
		}
		else
		{
		    if (end() || error())
		    {
			Report_Syntax_Error(Token(c, token_location),
					    "unexpected end of file; "
					    "did you forget a closing quote?");
		    }
		    else
		    {
			Check(c=='\n');
			Report_Syntax_Error(Token(c, token_location),
					    "unexpected end of line; "
					    "did you forget a closing quote?");
		    }
		}
	    }
	    text += '"';
	    pop_char();
	    Ensure(check_class_invariants());
	    return Token(STRING, text, token_location);
	}
	else   
	    // OTHER
	{
	    pop_char();
	    Ensure(check_class_invariants());
	    return Token(c, token_location);
	}

        Ensure(check_class_invariants());
        return Token(ERROR, token_location);
    }
}

//-----------------------------------------------------------------------//
/*! 
 * This function searches for the argument character in its internal
 * list of whitespace characters.  
 *
 * \param c
 * Character to be checked against the whitespace list.
 * 
 * \return \c true if and only if the character is found in the internal
 * whitespace list.
 */

bool Text_Token_Stream::is_whitespace(char const c) const
{
    return isspace(c) || whitespace_.count(c);
}

//-----------------------------------------------------------------------//
/*! 
 * This function searches for the argument character in the Token_Stream's
 * internal list of nonbreaking whitespace characters.  
 *
 * \param c
 * Character to be checked against the nonbreaking whitespace list.
 * 
 * \return \c true if and only if the character is found in the internal
 * nonbreaking whitespace list, and is \e not found in the breaking
 * whitespace list..
 */

bool Text_Token_Stream::is_nb_whitespace(char const c) const
{
    return !whitespace_.count(c) && (c==' ' || c=='\t');
}

//-------------------------------------------------------------------------//
/*!
 * An internal buffer is used to implement unlimited lookahead, necessary
 * for scanning numbers (which have a quite complex regular expression.)
 * This function pops a character off the top of the internal buffer, using
 * fill_character_buffer() if necessary to ensure that there is at least one
 * character in the buffer.   If the next character is a carriage return, the
 * line count is incremented.  
 *
 * \post <code>Result=='\n' && line==old_line+1  ||  line==old_line</code>
 *
 * \return The next character in the buffer.
 */

char Text_Token_Stream::pop_char()
{
    Remember(unsigned const old_line = line_);

    char const Result = peek();
    buffer_.pop_front();
    if (Result == '\n') line_++;

    Ensure(check_class_invariants());
    Ensure(Result == '\n' && line_ == old_line+1  ||  line_ == old_line);
    return Result;
}

//-------------------------------------------------------------------------//
/*!
 * \brief Try to scan a floating literal. 
 *
 * \return length of scanned literal; 0 if no literal could be scanned.
 */

unsigned Text_Token_Stream::scan_floating_literal()
{
    unsigned pos = 0;
    if (scan_fractional_constant(pos))
    {
	scan_exponent_part(pos);
	return pos;
    }
    else if (scan_digit_sequence(pos))
    {
	if (!scan_exponent_part(pos)) return 0;
	return pos;
    }
    else
    {
	return 0;
    }
}

//-------------------------------------------------------------------------//
/*!
 * \brief Try to scan a digit sequence. 
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_digit_sequence(unsigned &pos)
{
    unsigned const old_pos = pos;
    while (isdigit(peek(pos))) pos++;
    return pos-old_pos;
}   

//-------------------------------------------------------------------------//
/*!
 * \brief Try to scan an exponent part. 
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_exponent_part(unsigned &pos)
{
    unsigned const old_pos = pos;
    char c = peek(pos);
    if (c=='e' || c=='E')
    {
	pos++;
	c = peek(pos);
	if (c=='-' || c=='+') pos++;
	if (!scan_digit_sequence(pos))
	{
	    pos = old_pos;
	}
    }
    return pos-old_pos;
}   

//-------------------------------------------------------------------------//
/*!
 * \brief Try to scan a fractional constant. 
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_fractional_constant(unsigned &pos)
{
    unsigned const old_pos = pos;
    if (scan_digit_sequence(pos))
    {
	if (peek(pos) != '.')
	{
	    pos = old_pos;
	}
	else
	{
	    pos++;
	    scan_digit_sequence(pos);
	}

    }
    else if (peek(pos) == '.')
    {
	pos++;
	if (scan_digit_sequence(pos)==0)
	{
	    pos = old_pos;
	}
    }
    return pos-old_pos;
}

//-------------------------------------------------------------------------//
/*!
 * \brief Try to scan an integer literal. 
 *
 * \return length of scanned literal; 0 if no literal could be scanned.
 */

unsigned Text_Token_Stream::scan_integer_literal()
{
    unsigned pos = 0;
    if (scan_decimal_literal(pos))
    {
    }
    else if (scan_hexadecimal_literal(pos))
    {
    }
    else if (scan_octal_literal(pos))
    {
    }
    else
    {
	return 0;
    }
    return pos;
}

//-------------------------------------------------------------------------//
/*!
 * \brief Try to scan decimal literal. 
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_decimal_literal(unsigned &pos)
{
    unsigned const old_pos = pos;
    char c = peek(pos);
    if (isdigit(c) && c!='0')
    {
	while (isdigit(c))
	{
	    pos++;
	    c = peek(pos);
	}
    }
    return pos-old_pos;
}   

//-------------------------------------------------------------------------//
/*!
 * \brief Try to scan hexadecimal literal. 
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_hexadecimal_literal(unsigned &pos)
{
    unsigned old_pos = pos;
    if (peek(pos)=='0')
    {
	pos++;
	char c = peek(pos);
	if (c=='x' || c=='X')
	{
	    pos++;
	    c = peek(pos);
	    if (!isxdigit(c))
	    {
		pos = old_pos;
	    }
	    else
	    {
		while (isxdigit(c))
		{
		    pos++;
		    c = peek(pos);
		}
	    }
	}
	else
	{
	    pos = old_pos;
	}
    }
    return pos-old_pos;
}   

//-------------------------------------------------------------------------//
/*!
 * \brief Try to scan octal literal. 
 *
 * \return length of scanned text; 0 if no text could be scanned.
 */

unsigned Text_Token_Stream::scan_octal_literal(unsigned &pos)
{
    unsigned old_pos = pos;
    char c = peek(pos);
    if (c=='0')
    {
	while (isdigit(c) && c!='8' && c!='9')
	{
	    pos++;
	    c = peek(pos);
	}
    }
    return pos-old_pos;
}   

//-------------------------------------------------------------------------//
/*!
 * An internal buffer is used to implement unlimited lookahead,
 * necessary for scanning numbers (which have a quite complex regular
 * expression.)  This function peeks ahead the specified number of places
 * in the buffer, using fill_character_buffer() if necessary to ensure that
 * there is a sufficient number of characters in the buffer.
 *
 * \param pos
 * Position at which to peek.
 *
 * \return The character at \c buffer[pos].
 */

char Text_Token_Stream::peek(unsigned const pos)
{
    while (buffer_.size()<=pos)
    {
	fill_character_buffer();
    }

    Ensure(check_class_invariants());
    return buffer_[pos];
 }
  
//-------------------------------------------------------------------------//
/*!
 * This function flushes the Text_Token_Stream's internal buffers, so that
 * scanning resumes at the beginning of the file stream.  It is normally
 * called by its overriding version in children of Text_Token_Stream.
 */

void Text_Token_Stream::Rewind()
{
    buffer_.clear();
    line_ = 1;

    Token_Stream::Rewind();
    
    Ensure(check_class_invariants());
    Ensure(Error_Count()==0);
}

//---------------------------------------------------------------------------//
bool Text_Token_Stream::check_class_invariants() const
{
    return line_ > 0;
}

//---------------------------------------------------------------------------//
/* private */
void Text_Token_Stream::eat_whitespace()
{
    for (;;)
    {
	// Scan whitespace
	char c = peek();
	while (is_whitespace(c) && c!='\0')
	{
	    pop_char();
	    c = peek();
	}
	
	// Check for a comment
	if (c=='/')
	{
	    if (peek(1)=='/')
	    {
		// C++ comment
		while (c!='\n' && !error() && !end())
		{
		    pop_char();
		    c = peek();
		}
	    }
	    else if (peek(1)=='*')
	    {
		pop_char();  // pop the '/'
		pop_char();  // pop the '*'
		while ((peek(0) != '*' || peek(1) != '/') && !error() && 
		       !end())
		{
		    pop_char();
		}
		pop_char(); // pop the '*'
		pop_char(); // pop the '/'
	    }
	    else
	    {
		break;
	    }
	}
	else
	{
	    break;
	}
    }
    // private member function -- no invariant check
}
    
//---------------------------------------------------------------------//
/*!
 * \func Text_Token_Stream::location 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 *
 * This function returns a location string whose exact format is 
 * stream-specific.  For example, for a token stream that scans tokens
 * from a text file, this could be a string of the form "filename, line #".
 *
 * \return A string describing the location from which the 
 * Text_Token_Stream is currently scanning tokens.  
 */

//---------------------------------------------------------------------------//
/*! 
 * \param c Character to be pushed onto the back of the character queue.
 */

void Text_Token_Stream::character_push_back(char const c)
{
    buffer_.push_back(c);

    Ensure(check_class_invariants());
}

} // rtt_parser
//--------------------------------------------------------------------//
//                      end of Text_Token_Stream.cc
//--------------------------------------------------------------------//
