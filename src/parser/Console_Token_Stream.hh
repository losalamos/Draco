//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file Console_Token_Stream.hh
 * \author Kent G. Budge
 * \date Wed Jan 22 15:18:23 MST 2003
 * \brief Definition of class Console_Token_Stream.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef CCS4_Console_Token_Stream_HH
#define CCS4_Console_Token_Stream_HH

#include <fstream>
#include "Text_Token_Stream.hh"

namespace rtt_parser 
{
//-------------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Console-based token stream
 *
 * Console_Token_Stream represents a text token stream that derives its text
 * stream from the standard console input stream \c cin.  It reports errors
 * to the standard  console error stream \c cerr.
 *
 * This stream also differs from other streams in that the endline character
 * is converted to the semicolon character.  Parsers for use with console
 * streams typically treat the semicolon as an "end of statement" character
 * by specifying that it is NOT a whitespace character and looking for it as
 * a statement terminator.
 *
 * \note This class is an experimental concept and <i> should not be used in
 * production codes </i>.  In particular, the class cannot readily be tested
 * under our current unit testing system, since it is inherently interactive.
 *
 * \invariant <code>line > 0</code>
 */

class Console_Token_Stream : public Text_Token_Stream
{
  public:
    Console_Token_Stream();
    Console_Token_Stream(std::set<char> const &whitespace);
    
    void Rewind();
    
  protected:
    
    virtual std::string location() const;
    
    virtual void fill_character_buffer();
    virtual bool error() const;
    virtual bool end() const;
    
    virtual void Report_Error(const Token & token,
			      const std::string &message);
    
    virtual void Report_Error(const std::string &message);
};

} // rtt_parser

#endif  // CCS4_Console_Token_Stream_HH
//---------------------------------------------------------------------------//
//                      end of Console_Token_Stream.hh
//---------------------------------------------------------------------------//
