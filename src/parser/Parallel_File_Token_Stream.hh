//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file Parallel_File_Token_Stream.hh
 * \author Kent G. Budge
 * \date Wed Jan 22 15:18:23 MST 2003
 * \brief Definition of class Parallel_File_Token_Stream.
 * \note   Copyright © 2007 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef CCS4_Parallel_File_Token_Stream_HH
#define CCS4_Parallel_File_Token_Stream_HH

#include <fstream>
#include "Text_Token_Stream.hh"

namespace rtt_parser 
{
//-------------------------------------------------------------------------//
/*! 
 * \author Kent G. Budge
 * \date Thu Jan 23 08:41:54 MST 2003
 * \brief Parallel file-based token stream
 *
 * Parallel_File_Token_Stream is similar to File_Token_Stream.  However, it
 * assumes an SPMD (Single Program, Multiple Data) run environment.  Only the
 * designated I/O processor (normally processor 0) actually reads the
 * file. The characters read are then broadcast to the other processors.  The
 * advantage of parallelism at this level is that it avoids the I/O cost of
 * many processors reading one file while communicating data that is still
 * very flat.
 */

class Parallel_File_Token_Stream : public Text_Token_Stream
{
  public:

    //! Construct a Parallel_File_Token_Stream from a file.
    Parallel_File_Token_Stream(std::string const &filename);

    //! Construct a Parallel_File_Token_Stream from a file.
    Parallel_File_Token_Stream(std::string const &filename,
			       std::set<char> const &whitespace);
    
    void Rewind();
    
    virtual void Report(Token const & token,
                        std::string const &message);
    
    virtual void Report(std::string const &message);

    bool check_class_invariants() const;
    
  protected:
    
    virtual std::string location() const;
    
    //! Fill the character buffer.
    virtual void fill_character_buffer();

    virtual bool error() const;
    virtual bool end() const;

  private:

    // IMPLEMENTATION

    //! Open the input stream.
    void open();

    // DATA

    std::string filename_;  //!< File from which to take token text.
    std::ifstream infile_;  //!< Stream from which to take token text.

    bool is_io_processor_;     //!< Is this the designated I/O processor?

    bool at_eof_;        //!< Did processor 0 see the end of file?
    bool at_error_;      //!< Did processor 0 see an I/O error?
};

} // rtt_parser

#endif  // CCS4_Parallel_File_Token_Stream_HH
//---------------------------------------------------------------------------//
//                      end of Parallel_File_Token_Stream.hh
//---------------------------------------------------------------------------//
