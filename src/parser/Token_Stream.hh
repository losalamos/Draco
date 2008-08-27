//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file Token_Stream.hh
 * \author Kent G. Budge
 * \brief Definition of class Token_Stream.
 * \note Copyright 2006-2007 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef CCS4_Token_Stream_HH
#define CCS4_Token_Stream_HH

#include <deque>
#include <stdexcept>
#include "Token.hh"

namespace rtt_parser 
{
//-------------------------------------------------------------------------//
/*! 
 * \brief Parser exception class
 *
 * This is an exception class for reporting syntax errors in simple parsers.
 */

class Syntax_Error : public std::runtime_error
{
  public:

    // CREATORS
    
    Syntax_Error();
};

//-------------------------------------------------------------------------//
/*! 
 * \brief Abstract token stream for simple parsers
 *
 * A Token_Stream is a source of tokens for a Parse_Table or other parsing
 * code. Each token is returned as a Token struct.  There is unlimited
 * lookahead and pushback capability, but no backtracking is implemented in
 * this release.  In other words, the client may look as many tokens to the
 * right of the cursor as he desires, and he may push any number of tokens
 * onto the stream at the cursor position; but when the cursor advances (via
 * the Shift method) the token it advances over is discarded forever.
 */

class Token_Stream
    : private std::deque<Token>
{
  public:

    // CREATORS
    
    virtual ~Token_Stream(){}
    

    // MANIPULATORS

    //! Return the next token in the stream and advance the cursor.
    Token shift();
    
    //! Look ahead at tokens.
    Token lookahead(unsigned pos=0);
    
    //! Insert a token into the stream at the cursor position.
    void pushback(Token const &token);

    //-----------------------------------------------------------------------//
    /*! 
     * \brief Reset the stream
     *
     * This function resets the token stream to some initial state defined
     * by the child class.
     */
    virtual void rewind() = 0;

    virtual void report_syntax_error(Token const &token,
				     std::string const &message);

    virtual void report_syntax_error(std::string const &message);
    
    virtual void report_semantic_error(Token const &token,
				       std::string const &message);
    
    virtual void report_semantic_error(std::string const &message);
    virtual void report_semantic_error(std::exception const &message);

    
    //-----------------------------------------------------------------------//
    /*! 
     * \author Kent G. Budge
     * \date Thu Jan 23 08:41:54 MST 2003
     * \brief Report an error to the user.
     *
     * This function sends a message to the user in a stream-specific manner.
     *
     * \param token
     * Token with which the message is associated.
     * \param message
     * Message to be passed to the user.
     */
    virtual void report(Token const &token,
                                    std::string const &message) = 0;
    
    //-----------------------------------------------------------------------//
    /*! 
     * \author Kent G. Budge
     * \date Thu Jan 23 08:41:54 MST 2003
     * \brief Report an error to the user.
     *
     * This function sends a message to the user in a stream-specific
     * manner.  This variant assumes that the cursor gives the correct
     * location.
     *
     * \param message
     * Message to be passed to the user.
     */
    virtual void report(std::string const &message) = 0;
    

    // ACCESSORS

    //! Return the number of errors reported to the stream. 
    unsigned error_count() const { return error_count_; }

    //! The current implementation of Token_Stream has no invariants.
    bool check_class_invariants() const { return true; }
    
    // STATICS
    
  protected:

    // IMPLEMENTATION

    //! Construct a Token_Stream.
    inline Token_Stream();

    //-----------------------------------------------------------------------//
    /*! 
     * \author Kent G. Budge
     * \date Thu Jan 23 08:41:54 MST 2003
     * \brief Add one or more tokens to the end of the lookahead buffer.
     *
     * This function must be overridden by all child classes to scan 
     * tokens from the ultimate token source (such as a text file or GUI).
     * If no more tokens are available, the overriding function must return
     * an EXIT token onto the end of the token buffer.
     *
     * \return Next token to put on the token buffer.
     */
    
    virtual Token fill_() = 0;
    
  private:

    // DATA

    unsigned error_count_;  //!< Number of errors reported.
};

//-----------------------------------------------------------------------//
/*! 
 * Construct a Token_Stream and place the cursor at the start of the
 * stream.
 */

inline Token_Stream::Token_Stream() 
    : error_count_(0) 
{
    Ensure(check_class_invariants());
    Ensure(error_count() == 0);
}

}  // namespace rtt_parser

#endif  // CCS4_Token_Stream_HH
//---------------------------------------------------------------------------//
//                          end of Token_Stream.hh
//---------------------------------------------------------------------------//
