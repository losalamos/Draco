//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/Ensight_Stream.cc
 * \author Rob Lowrie
 * \date   Mon Nov 15 10:03:51 2004
 * \brief  Ensight_Stream implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Ensight_Stream.hh"

#include <iomanip>
#include <ds++/Assert.hh>
#include <ds++/Packing_Utils.hh>

namespace rtt_viz
{

//---------------------------------------------------------------------------//
/*!
 * \brief The endl manipulator.
 *
 * Note that this is a function within the rtt_viz namespace, NOT a member
 * function of Ensight_Stream.
 */
Ensight_Stream& endl(Ensight_Stream &s)
{
    if ( ! s.d_binary )
	s.d_stream << std::endl;
    
    return s;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * \a geom_file is used only so that the "C Binary" header may be dumped when
 * \a binary is true.  If the geometry file is binary, Ensight assumes that
 * all data files are also binary.  This class does NOT check whether \a
 * binary is consistent across all geometry and data files.
 *
 * \param file_name  Name of output file.   
 * \param binary     If true, output binary.  Otherwise, output ascii.
 * \param geom_file  If true, then a geometry file will be dumped.
 */
Ensight_Stream::Ensight_Stream(const std::string &file_name,
			       const bool binary,
			       const bool geom_file)
    : d_binary(binary)
{
    Require(! file_name.empty());
    
    // Open the stream.
    
    if ( binary )
	d_stream.open(file_name.c_str(), std::ios::binary);
    else
	d_stream.open(file_name.c_str());

    Check(d_stream);

    // Set up the file.
	
    if ( binary )
    {
	if ( geom_file )
	    *this << "C Binary";
    }
    else
    {
	// set precision for ascii mode
	d_stream.precision(5);
	d_stream.setf(std::ios::scientific, std::ios::floatfield);
    }

    Ensure(d_stream.good());
}

//---------------------------------------------------------------------------//
/*!
 * \brief Output for ints.
 */
Ensight_Stream& Ensight_Stream::operator<<(const int i)
{
    if ( d_binary )
	binary_write(i);
    else
	d_stream << std::setw(10) << i;
    
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Output for size_t.
 *
 * This is a convience function.  It simply casts to int.  Ensight does not
 * support output of unsigned ints.
 */
Ensight_Stream& Ensight_Stream::operator<<(const std::size_t i)
{
    int j(i);
    Check(j >= 0);
    *this << j;
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Output for doubles.
 *
 * Note that Ensight only supports "float" for binary mode.
 */
Ensight_Stream& Ensight_Stream::operator<<(const double d)
{
    if ( d_binary )
	binary_write(float(d));
    else
	d_stream << std::setw(12) << d;
    
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Output for strings.
 */
Ensight_Stream& Ensight_Stream::operator<<(const std::string &s)
{
    if ( d_binary )
    {
	// Ensight demands all character strings be 80 chars.  Make it so.
	std::string sc(s);
	sc.resize(80);
	d_stream.write(sc.c_str(), 80);
    }
    else
	d_stream << s;
    
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Output for function pointers.
 */
Ensight_Stream& Ensight_Stream::operator<<(FP f)
{
    Require(f);
    
    f(*this);
    return *this;
}

//---------------------------------------------------------------------------//
// PRIVATE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Does binary write of \a v.
 *
 * The type \a T must support sizeof(T).
 */
// The template implementation is defined here because only functions within
// this translation unit should be calling this function.
template <class T>
void Ensight_Stream::binary_write(const T v)
{
    char *vc = new char[sizeof(T)];

    rtt_dsxx::Packer p;
    p.set_buffer(sizeof(T), vc);
    p.pack(v);

    d_stream.write(vc, sizeof(T));
    delete[] vc;
}

} // end of rtt_viz

//---------------------------------------------------------------------------//
//                              end of Ensight_Stream.cc
//---------------------------------------------------------------------------//
