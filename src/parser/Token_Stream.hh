//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file Token_Stream.hh
 * \author Kent G. Budge
 * \date Wed Jan 22 13:15:29 MST 2003
 * \brief Definition of class Token_Stream.
 * \note Copyright @ 2003 The Regents of the University of California.
 *
 * revision history:
 * 0) original
 * 1) kgbudge (03/12/03): Remove erroneous precondition from Report_Error.
 * 2) kgbudge (03/08/10): Solo inspection of documentation, assertions, and
 * tests. 
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
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Parser exception class
 *
 * This is an exception class for reporting syntax errors in simple parsers.
 */

class Syntax_Error : public std::runtime_error
{
  public:
    Syntax_Error();
};

//-------------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Abstract token stream for simple parsers
 *
 * A Token_Stream is a source of tokens for a Parse_Table or other parsing
 * code. Each token is returned as a Token struct.  There is unlimited
 * lookahead and pushback capability, but no backtracking is implemented in
 * this release.  In other words, the client may look as many tokens to the
 * right of the cursor as he desires, and he may push any number of tokens
 * onto the stream at the cursor position; but when the cursor advances (via
 * the Shift method) the token it advances over is discarded.
 */

class Token_Stream : private std::deque<Token>
{
  protected:

    //! Construct a Token_Stream.
    inline Token_Stream();
    
  public:

    //! This is an abstract class.
    virtual ~Token_Stream(){}
    
    Token Shift();
    
    //! Look ahead at tokens.
    Token Lookahead(unsigned pos=0);
    
    //! Insert a token into the stream at the cursor position.
    void Pushback(Token const &token);
    
    virtual void Report_Syntax_Error(Token const &token,
				     std::string const &message);

    virtual void Report_Syntax_Error(std::string const &message);
    
    virtual void Report_Semantic_Error(Token const &token,
				       std::string const &message);
    
    virtual void Report_Semantic_Error(std::string const &message);
    virtual void Report_Semantic_Error(std::exception const &message);
    
    //! Return the number of errors reported to the stream. 
    unsigned Error_Count() const { return error_count; }

    //! The current implementation of Token_Stream has no invariants.
    bool check_class_invariants() const { return true; }

    //! The current implementation of Token_Stream has no static invariants.
    bool check_static_class_invariants() const { return true; }
    
  protected:
    
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
    
    virtual Token fill() = 0;
    
    //-----------------------------------------------------------------------//
    /*! 
     * \author Kent G. Budge
     * \date Thu Jan 23 08:41:54 MST 2003
     * \brief Report an error to the user.
     *
     * This function reports an error to the user in a stream-specific
     * manner.
     *
     * \param token
     * Token at which the error occurred.
     * \param message
     * Message to be passed to the user.
     */
    virtual void Report_Error(Token const &token,
			      std::string const &message) = 0;
    
    //-----------------------------------------------------------------------//
    /*! 
     * \author Kent G. Budge
     * \date Thu Jan 23 08:41:54 MST 2003
     * \brief Report an error to the user.
     *
     * This function reports an error to the user in a stream-specific
     * manner.  This variant assumes that the cursor gives the correct
     * location.
     *
     * \param message
     * Message to be passed to the user.
     */
    virtual void Report_Error(std::string const &message) = 0;

//---------------------------------------------------------------------------//
/*! 
 * \brief Reset the stream
 *
 * This function resets the token stream to some initial state defined
 * by the child class.
 *
 * \post \c Error_Count()==0
 */

    virtual void Rewind() = 0;
    
  private:

    unsigned error_count;  //!< Number of errors reported.

    Token_Stream(Token_Stream const &);            //!< Not copyable
    Token_Stream& operator=(Token_Stream const &); //!< Not assignable
};

//-----------------------------------------------------------------------//
/*! 
 * Construct a Token_Stream and place the cursor at the start of the
 * stream.
 *
 * \post <code>Error_Count() == 0</code>
 */

inline Token_Stream::Token_Stream() 
    : error_count(0) 
{
    Ensure(check_class_invariants());
    Ensure(Error_Count() == 0);
}

}  // namespace rtt_parser

#endif  // CCS4_Token_Stream_HH
//---------------------------------------------------------------------------//
//                          end of Token_Stream.hh
//---------------------------------------------------------------------------//
