//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file File_Token_Stream.hh
 * \author Kent G. Budge
 * \date Wed Jan 22 15:18:23 MST 2003
 * \brief Definition of class File_Token_Stream.
 * \note   Copyright © 2003 The Regents of the University of California.
 *
 * revision history:
 * 0) original
 * 1) kgbudge (03/08/10): 
 *    Solo inspection of documentation, assertions, and tests. 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef CCS4_File_Token_Stream_HH
#define CCS4_File_Token_Stream_HH

#include <fstream>
#include "Text_Token_Stream.hh"

namespace rtt_parser 
{
//-------------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief File-based token stream
 *
 * File_Token_Stream represents a text token stream that derives its text
 * stream from a file in the filesystem.  It reports errors to the standard 
 * console error stream \c cerr.
 *
 * \invariant <code>line > 0</code>
 */

class File_Token_Stream : public Text_Token_Stream
{
  public:
    File_Token_Stream(std::string filename);
    File_Token_Stream(std::string filename, std::set<char> const &whitespace);

    void open(std::string filename);
    
    void Rewind();
    
  protected:
    
    virtual std::string location() const;
    
    virtual void fill_character_buffer();
    virtual bool error() const;
    virtual bool end() const;
    
    virtual void Report_Error(Token const & token,
			      std::string const &message);
    
    virtual void Report_Error(std::string const &message);

  private:
    std::string filename;  //!< File from which to take token text.
    std::ifstream infile;  //!< Stream from which to take token text.
};

} // rtt_parser

#endif  // CCS4_File_Token_Stream_HH
//---------------------------------------------------------------------------//
//                      end of File_Token_Stream.hh
//---------------------------------------------------------------------------//
