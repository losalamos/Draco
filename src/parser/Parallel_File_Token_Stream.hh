//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   Parallel_File_Token_Stream.hh
 * \author Kent G. Budge
 * \brief  Definition of class Parallel_File_Token_Stream.
 * \note   Copyright (C) 2006-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef CCS4_Parallel_File_Token_Stream_HH
#define CCS4_Parallel_File_Token_Stream_HH

#include <fstream>
#include "Text_Token_Stream.hh"

#if defined(MSVC)
#   pragma warning (push)
#   pragma warning (disable:4251) // warning C4251: 'rtt_parser::Parallel_File_Token_Stream::filename_' : class 'std::basic_string<_Elem,_Traits,_Ax>' needs to have dll-interface to be used by clients of class 'rtt_parser::Parallel_File_Token_Stream'
#endif

namespace rtt_parser 
{
using std::string;
using std::set;
using std::ifstream;

//-------------------------------------------------------------------------//
/*! 
 * \brief Parallel file-based token stream
 *
 * \c Parallel_File_Token_Stream is similar to \c File_Token_Stream. However,
 * it assumes an SPMD (Single Program, Multiple Data) run environment. Only
 * the designated I/O processor (normally processor 0) actually reads the
 * file. The characters read are then broadcast to the other processors. The
 * advantage of parallelism at this level is that it avoids the I/O cost of
 * many processors reading one file while communicating data that is still
 * very flat.
 */

class DLL_PUBLIC Parallel_File_Token_Stream : public Text_Token_Stream
{
  public:

    // CREATORS

    //! Construct a Parallel_File_Token_Stream from a file.
    Parallel_File_Token_Stream(string const &filename);

    //! Construct a Parallel_File_Token_Stream from a file.
    Parallel_File_Token_Stream(string const &filename,
			       set<char> const &whitespace);

    // MANIPULATORS

    //! Rewind the Parallel_File_Token_Stream.
    virtual void rewind();

    //! Report a condition.
    virtual void report(Token const & token,
                        string const &message);

    //! Report a condition.
    virtual void report(string const &message);

    // ACCESSORS

    //! Check the class invariants.
    bool check_class_invariants() const;
    
  protected:

    // IMPLEMENTATION
    
    virtual string location_() const;
    
    virtual void fill_character_buffer_();

    virtual bool error_() const;
    virtual bool end_() const;

  private:

    //! Open the input stream.
    void open_();

    // DATA

    string filename_;  //!< File from which to take token text.
    ifstream infile_;  //!< Stream from which to take token text.

    bool is_io_processor_;     //!< Is this the designated I/O processor?

    bool at_eof_;        //!< Did processor 0 see the end of file?
    bool at_error_;      //!< Did processor 0 see an I/O error?
};

} // rtt_parser

#if defined(MSVC)
#   pragma warning (pop)
#endif

#endif  // CCS4_Parallel_File_Token_Stream_HH
//---------------------------------------------------------------------------//
// end of Parallel_File_Token_Stream.hh
//---------------------------------------------------------------------------//
