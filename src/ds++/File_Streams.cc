//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/File_Streams.cc
 * \author Rob Lowrie
 * \date   Mon Nov 15 10:03:51 2004
 * \brief  File_Streams implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "File_Streams.hh"
#include <iomanip>

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
// File_Output functions.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param filename The file name to open for writing.  If empty, open() must
 *                 be used later to open a file.
 * \param binary   If true, use binary mode for writing.
 */
File_Output::
File_Output(const std::string &filename,
	    const bool binary)
    : d_last_was_char(false)
    , d_binary(binary)
  
{
    if ( ! filename.empty() ) open(filename, binary);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor
 */
File_Output::~File_Output()
{
    close();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opens a file for writing.
 *
 * \param filename The file name to open for writing.
 * \param binary   If true, use binary mode for writing.
 */
void File_Output::open(const std::string &filename,
		       const bool binary)
{
    Require(! filename.empty());

    if ( d_stream.is_open() ) close();
    
    d_last_was_char = false;
    d_binary = binary;

    if ( d_binary )
	d_stream.open(filename.c_str(), std::ios::binary);
    else
	d_stream.open(filename.c_str());

    Ensure(d_stream);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Closes the stream.
 */
void File_Output::close()
{
    if ( d_stream.is_open() )
    {
	if ( d_last_was_char ) d_stream << std::endl;
	d_stream.close();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Stream output for type char.
 */
File_Output& File_Output::operator<<(const char c)
{
    Require(d_stream.is_open());
    
    if ( d_binary )
    {
	d_stream.write(&c, 1);
    }
    else // ascii mode
    {
	// For char, we don't add a newline, in case its part of a
	// character string.
	d_last_was_char = true;
	d_stream << c;
    }

    Ensure(d_stream.good());
    
    return *this;
}

//---------------------------------------------------------------------------//
// File_Input functions.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \param filename The file name to open for reading.  If empty, open() must
 *                 be used later to open a file.
 * \param binary   If true, use binary mode for reading.
 */
File_Input::
File_Input(const std::string &filename,
	   const bool binary)
    : d_char_line(-1)
    , d_binary(binary)
{
    if ( ! filename.empty() ) open(filename, binary);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor
 */
File_Input::~File_Input()
{
    close();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opens a file for reading.
 *
 * \param filename The file name to open for reading.
 * \param binary   If true, use binary mode for reading.
 */
void File_Input::open(const std::string &filename,
		      const bool binary)
{
    Require(! filename.empty());

    d_binary = binary;
    d_char_line = -1;

    if ( d_stream.is_open() ) close();

    if ( d_binary )
    {
	d_stream.open(filename.c_str(), std::ios::binary);
    }
    else
    {
	d_stream.open(filename.c_str());
    }

    Ensure(d_stream);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Closes the stream.
 */
void File_Input::close()
{
    if ( d_stream.is_open() )
    {
	d_stream.close();
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Stream input for type char.
 */
File_Input& File_Input::operator>>(char &c)
{
    Require(d_stream.is_open());
    
    if ( d_binary )
    {
	d_stream.read(&c, 1);
    }
    else // ascii mode
    {
	if ( d_char_line < 0 )
	{
	    std::getline(d_stream, d_line);
	    Check(! d_line.empty());
	    d_char_line = 0;
	}

	Check(d_char_line < d_line.size());
	c = d_line[d_char_line];
	++d_char_line;
    }

    Ensure(d_stream.good());
    
    return *this;
}

} // end of rtt_dsxx

//---------------------------------------------------------------------------//
//                              end of ds++/File_Streams.cc
//---------------------------------------------------------------------------//
