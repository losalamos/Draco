//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/Ensight_Stream.cc
 * \author Rob Lowrie
 * \date   Mon Nov 15 10:03:51 2004
 * \brief  Ensight_Stream implementation file.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Ensight_Stream.hh"
#include "ds++/Assert.hh"
#include "ds++/Packing_Utils.hh"
#include <iomanip>

namespace rtt_viz {

//---------------------------------------------------------------------------//
/*!
 * \brief The endl manipulator.
 *
 * Note that this is a function within the rtt_viz namespace, NOT a member
 * function of Ensight_Stream.
 */
Ensight_Stream &endl(Ensight_Stream &s) {
  Require(s.d_stream.is_open());

  if (!s.d_binary)
    s.d_stream << '\n';

  Require(s.d_stream.good());

  return s;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor
 *
 * This constructor opens the stream, if \a file_name is non-empty.  See open()
 * for more information.
 *
 * \param file_name  Name of output file.
 * \param binary     If true, output binary.  Otherwise, output ascii.
 * \param geom_file  If true, then a geometry file will be dumped.
 */
Ensight_Stream::Ensight_Stream(const std::string &file_name, const bool binary,
                               const bool geom_file)
    : d_stream(), d_binary(false) {
  if (!file_name.empty())
    open(file_name, binary, geom_file);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor
 *
 * Automatically closes stream, if open.
 */
Ensight_Stream::~Ensight_Stream(void) { close(); }

//---------------------------------------------------------------------------//
/*!
 * \brief Opens the stream.
 *
 * \a geom_file is used only so that the "C Binary" header may be dumped when
 * \a binary is true.  If the geometry file is binary, Ensight assumes that all
 *         data files are also binary.  This class does NOT check whether \a
 *         binary is consistent across all geometry and data files.
 *
 * \param file_name  Name of output file.
 * \param binary     If true, output binary.  Otherwise, output ascii.
 * \param geom_file  If true, then a geometry file will be dumped.
 */
void Ensight_Stream::open(const std::string &file_name, const bool binary,
                          const bool geom_file) {
  Require(!file_name.empty());

  d_binary = binary;

  // Open the stream.
  if (binary)
    d_stream.open(file_name.c_str(), std::ios::binary);
  else
    d_stream.open(file_name.c_str());

  Check(d_stream);

  // Set up the file.

  if (binary) {
    if (geom_file)
      *this << "C Binary";
  } else {
    // set precision for ascii mode
    d_stream.precision(5);
    d_stream.setf(std::ios::scientific, std::ios::floatfield);
  }

  Ensure(d_stream.good());
}

//---------------------------------------------------------------------------//
//! Closes the stream.
void Ensight_Stream::close() {
  if (d_stream.is_open()) {
    d_stream.close();
  }
}

//---------------------------------------------------------------------------//
//! Output for ints.
Ensight_Stream &Ensight_Stream::operator<<(const int32_t i) {
  Require(d_stream.is_open());

  if (d_binary)
    binary_write(i);
  else
    d_stream << std::setw(10) << i;

  Ensure(d_stream.good());

  return *this;
}

//---------------------------------------------------------------------------//
/*
 * \brief Output for unsigned
 *
 * This is a convience function.  It simply casts to int.  Ensight does not
 * support output of unsigned ints.
 */
Ensight_Stream &Ensight_Stream::operator<<(const unsigned i) {
  Check(i < INT_MAX);
  int const j = static_cast<int>(i);
  *this << j;
  return *this;
}

//---------------------------------------------------------------------------//
/*
 * \brief Output for int64_t.
 *
 * This is a convience function.  It simply casts to int.  Ensight does not
 * support output of unsigned ints.
 *
 * \bug Not tested so commented out.
 */
// Ensight_Stream &Ensight_Stream::operator<<(const int64_t i) {
//   Check(i < INT_MAX && i > -1 * INT_MAX);
//   int const j = static_cast<int>(i);
//   *this << j;
//   return *this;
// }

//---------------------------------------------------------------------------//
/*!
 * \brief Output for uint64_t.
 *
 * This is a convience function.  It simply casts to int.  Ensight does not
 * support output of unsigned ints.
 *
 * \bug Not tested so commented out.
 */
// Ensight_Stream &Ensight_Stream::operator<<(const uint64_t i) {
//   Check(i < INT_MAX);
//   int const j = static_cast<int>(i);
//   *this << j;
//   return *this;
// }

//---------------------------------------------------------------------------//
/*!
 * \brief Output for doubles.
 *
 * Note that Ensight only supports "float" for binary mode.
 */
Ensight_Stream &Ensight_Stream::operator<<(const double d) {

#if defined(MSVC) && MSVC_VERSION < 1900
  // [2015-02-06 KT]: By default, MSVC uses a 3-digit exponent (presumably
  // because numeric_limits<double>::max() has a 3-digit exponent.)
  // Enable two-digit exponent format to stay consistent with GNU and
  // Intel on Linux.(requires <stdio.h>).
  unsigned old_exponent_format = _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  Require(d_stream.is_open());

  if (d_binary)
    binary_write(float(d));
  else
    d_stream << std::setw(12) << d;

  Ensure(d_stream.good());

#if defined(MSVC) && MSVC_VERSION < 1900
  // Disable two-digit exponent format
  _set_output_format(old_exponent_format);
#endif

  return *this;
}

//---------------------------------------------------------------------------//
//! Output for strings.
Ensight_Stream &Ensight_Stream::operator<<(const std::string &s) {
  Require(d_stream.is_open());

  if (d_binary) {
    // Ensight demands all character strings be 80 chars.  Make it so.
    std::string sc(s);
    sc.resize(80);
    d_stream.write(sc.c_str(), 80);
  } else
    d_stream << s;

  Ensure(d_stream.good());

  return *this;
}

//---------------------------------------------------------------------------//
//! Output for function pointers.
Ensight_Stream &Ensight_Stream::operator<<(FP f) {
  Require(d_stream.is_open());

  Require(f);

  f(*this);

  Ensure(d_stream.good());

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
 *
 * The template implementation is defined here because only functions within
 * this translation unit should be calling this function.
 */
template <typename T> void Ensight_Stream::binary_write(const T v) {
  Require(d_stream.is_open());

  char *vc = new char[sizeof(T)];

  rtt_dsxx::Packer p;
  p.set_buffer(sizeof(T), vc);
  p.pack(v);

  d_stream.write(vc, sizeof(T));
  delete[] vc;

  Ensure(d_stream.good());
}

} // namespace rtt_viz

//---------------------------------------------------------------------------//
// end of Ensight_Stream.cc
//---------------------------------------------------------------------------//
