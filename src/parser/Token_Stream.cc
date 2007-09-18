//----------------------------------*-C++-*----------------------------------//
/*!
 * \file Token_Stream.cc
 * \author Kent G. Budge
 * \date Wed Jan 22 14:57:49 MST 2003
 * \brief Definitions of Token_Stream member functions.
 * \note   Copyright © 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Token_Stream.hh"

namespace rtt_parser 
{
using namespace std;

//---------------------------------------------------------------------------//
/*! 
 * \brief Construct a Syntax_Error exception object.
 *
 * \post \c what()=="syntax error"
 */

Syntax_Error::Syntax_Error()
    : runtime_error("syntax error")
{
    using namespace std;

    Ensure(!strcmp(what(), "syntax error"));
}

//-----------------------------------------------------------------------//
/*! 
 * \brief Return the next token in the Token_Stream.
 *
 * This function returns the token at the cursor position and advance the 
 * cursor. It will, if necessary, fill() the token buffer first.
 *
 * \return <code>old Lookahead()</code>
 */

Token Token_Stream::Shift()
{    
    Token const Result = Lookahead();
    pop_front();
    
    Ensure(check_class_invariants());
    // Ensure the cursor advances one place to the right, discarding the
    // leftmost token.
    return Result;
}

//-----------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 *
 * This function looks ahead in the token stream without changing the
 * cursor position.  It will, if necessary, fill() the token buffer
 * first.  If the requested position is at or past the end of the file, an
 * EXIT token will be returned.
 *
 * \param pos Number of tokens to look ahead, with 0 being the token at the
 * cursor position.
 *
 * \return The token at the specified position relative to the
 * cursor.
 */

Token Token_Stream::Lookahead(unsigned const pos)
{
    while (size()<=pos)
    {
	push_back(fill());
    }

    Ensure(check_class_invariants());
    return operator[](pos);
}

//-----------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 *
 * This function pushes the specified token onto the front of the
 * token stream, so that it is now the token in the Lookahead(0)
 * position.
 */

void Token_Stream::Pushback(Token const &token)
{
    push_front(token);

    Ensure(check_class_invariants());
    Ensure(Lookahead()==token);
}

//-----------------------------------------------------------------------//
/*! 
 * \brief Report a syntax error to the user.
 *
 * The default implementation of this function passes its message
 * on to Report_Error, then throws a Syntax_Error exception.
 *
 * A syntax error is a badly formed construct that requires explicit
 * error recovery (including resynchronization) by the parsing software.
 *
 * \param token
 * Token at which the error occurred.
 * \param message
 * The message to be passed to the user.
 *
 * \throw Syntax_Error This function never returns.  It always throws a
 * Syntax_Error exception to be handled by the parsing software.
 */

void Token_Stream::Report_Syntax_Error(Token const &token,
				       string const &message)
{
    try 
    {
	error_count_++;
	Report(token, message);
    }
    catch(...)
    {
	// An error at this point really hoses us.  It means something went
	// sour with the reporting mechanism, and there probably isn't much
	// we can do about it.
	throw std::bad_exception();
    }
    
    Ensure(check_class_invariants());
    throw Syntax_Error();
}

//-----------------------------------------------------------------------//
/*! 
 * \brief Report a syntax error to the user.
 *
 * The default implementation of this function passes its message
 * on to Report, then throws a Syntax_Error exception.
 *
 * A syntax error is a badly formed construct that requires explicit
 * error recovery (including resynchronization) by the parsing software.
 *
 * This versiona ssumes that the cursor is the location of the error.
 *
 * \param message
 * The message to be passed to the user.
 *
 * \post <code>Error_Count()==old Error_Count()+1</code>
 *
 * \throw Syntax_Error This function never returns.  It always throws a
 * Syntax_Error exception to be handled by the parsing software.
 */

void Token_Stream::Report_Syntax_Error(string const &message)
{
    try 
    {
	error_count_++;
	Report(message);
    }
    catch(...)
    {
	// An error at this point really hoses us.  It means something went
	// sour with the reporting mechanism, and there probably isn't much
	// we can do about it.
	throw std::bad_exception();
    }
    
    throw Syntax_Error();
}

//-----------------------------------------------------------------------//
/*! 
 * \brief Report a semantic error to the user.
 *
 * The default implementation of this function passes its message
 * on to Report, then returns.
 *
 * A semantic error is a well-formed construct that has a bad value.  Because
 * the construct is well-formed, parsing may be able to continue after the
 * error is reported without any explicit recovery by the parsing software.
 *
 * \param token
 * Token at which the error occurred.
 * \param message
 * The message to be passed to the user.
 */

void Token_Stream::Report_Semantic_Error(Token const &token,
					 string const &message)
{
    error_count_++;
    Report(token, message);

    Ensure(check_class_invariants());
}

//-----------------------------------------------------------------------//
/*! 
 * \brief Report a semantic error to the user.
 *
 * The default implementation of this function passes its message
 * on to Report, then returns.
 *
 * A semantic error is a well-formed construct that has a bad value.  Because
 * the construct is well-formed, parsing may be able to continue after the
 * error is reported without any explicit recovery by the parsing software.
 *
 * This version assumes that the cursor is the error location.
 *
 * \param message
 * The message to be passed to the user.
 *
 * \post <code>Error_Count()==old Error_Count()+1</code>
 */

void Token_Stream::Report_Semantic_Error(string const &message)
{
    error_count_++;
    Report(message);

    Ensure(check_class_invariants());
}

//-----------------------------------------------------------------------//
/*! 
 * \brief Report a semantic error to the user.
 *
 * The default implementation of this function passes its message
 * on to Report, then returns.
 *
 * A semantic error is a well-formed construct that has a bad value.  Because
 * the construct is well-formed, parsing may be able to continue after the
 * error is reported without any explicit recovery by the parsing software.
 *
 * This version assumes that the cursor is the error location.
 *
 * \param message
 * The exception whose message is to be passed to the user.
 *
 * \post <code>Error_Count()==old Error_Count()+1</code>
 */

void Token_Stream::Report_Semantic_Error(exception const &message)
{
    error_count_++;
    Report(message.what());

    Ensure(check_class_invariants());
}
    
//---------------------------------------------------------------------------//
/*! 
 * \brief Reset the token stream.
 *
 * This function is normally called by its overriding version in children of
 * Token_Stream. It flushes the queues and resets the error count.
 */

void Token_Stream::Rewind()
{
    Require(check_class_invariants());

    error_count_ = 0;
    clear();

    Ensure(check_class_invariants());
}

} // rtt_parser
//---------------------------------------------------------------------------//
//                          end of Token_Stream.cc
//---------------------------------------------------------------------------//
