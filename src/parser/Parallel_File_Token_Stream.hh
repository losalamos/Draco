//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Parallel_File_Token_Stream.hh
 * \author Kent G. Budge
 * \brief  Definition of class Parallel_File_Token_Stream.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#ifndef CCS4_Parallel_File_Token_Stream_HH
#define CCS4_Parallel_File_Token_Stream_HH

#include "Text_Token_Stream.hh"
#include <fstream>
#include <vector>

namespace rtt_parser {
using std::ifstream;
using std::set;
using std::string;

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

class DLL_PUBLIC_parser Parallel_File_Token_Stream : public Text_Token_Stream {
public:
  // CREATORS

  //! Construct an empty Parallel_File_Token_Stream
  Parallel_File_Token_Stream();

  //! Construct a Parallel_File_Token_Stream from a file.
  Parallel_File_Token_Stream(string const &filename);

  //! Construct a Parallel_File_Token_Stream from a file.
  Parallel_File_Token_Stream(string const &filename,
                             set<char> const &whitespace);

  // MANIPULATORS

  //! Reopen the Parallel_File_Token_Stream with a new file.
  void open(string const &filename);

  //! Rewind the Parallel_File_Token_Stream.
  virtual void rewind();

  //! Report a condition.
  virtual void report(Token const &token, string const &message);

  //! Report a condition.
  virtual void report(string const &message);

  //! Report a comment.
  virtual void comment(std::string const &message);

  // ACCESSORS

  //! Check the class invariants.
  bool check_class_invariants() const;

protected:
  // IMPLEMENTATION

  virtual string location_() const;

  virtual void fill_character_buffer_();

  virtual bool error_() const;
  virtual bool end_() const;

  virtual void push_include(std::string &include_file_name);
  virtual void pop_include();

private:
  // NESTED TYPES

  struct letter {
    // IMPLEMENTATION

    //! Constructor
    letter(string const &file_name);

    bool check_class_invariants() const;

    //! Open the input stream.
    void open_();

    //! Read the next set of characters into a buffer.
    void fill_character_buffer_(std::vector<char> &comm_buffer);

    //! Rewind the stream.
    void rewind();

    // DATA
    string filename_; //!< File from which to take token text.
    ifstream infile_; //!< Stream from which to take token text.

    bool is_io_processor_; //!< Is this the designated I/O processor?

    bool at_eof_;   //!< Did processor 0 see the end of file?
    bool at_error_; //!< Did processor 0 see an I/O error?
  };

  // IMPLEMENTATION

  // DATA

  std::stack<std::shared_ptr<letter>> letters_;
  std::shared_ptr<letter> letter_;
};

} // namespace rtt_parser

#endif // CCS4_Parallel_File_Token_Stream_HH

//---------------------------------------------------------------------------//
// end of Parallel_File_Token_Stream.hh
//---------------------------------------------------------------------------//
