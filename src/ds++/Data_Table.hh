//----------------------------------*-C++-*----------------------------------//
/*!
 * \file    Data_Table.hh
 * \author  Paul Henning
 * \brief   Declaration of class Data_Table
 * \note    Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 * \version $Id$
 */
//---------------------------------------------------------------------------//
#ifndef dsxx_Data_Table_hh
#define dsxx_Data_Table_hh

#include "Assert.hh"
#include "ds++/config.h"
#include <vector>

/*!
  Data_Table provides read-only, DBC-checked, container-like
  access to a sequential range of memory locations, or a scalar.  This is
  useful in situations where the amount of data changes depending on
  compile-time factors, but you always want to access it as an array.
*/

namespace rtt_dsxx {

//! Provide const array-style access to an actual array or a scalar.
template <typename T> class Data_Table {
public:
  typedef T const *const_iterator;

public:
  inline Data_Table(Data_Table const &);
  inline explicit Data_Table(std::vector<T> const &v);
  inline Data_Table(const_iterator const begin, const_iterator const end);
  inline explicit Data_Table(T const &value);
  inline Data_Table();
  inline T const &operator[](const unsigned i) const;
  inline const_iterator begin() const { return d_begin; }
  inline const_iterator end() const { return d_end; }
  inline unsigned size() const { return d_end - d_begin; }
  inline T const &front() const;
  inline T const &back() const;
  inline T *access();
  Data_Table &operator=(Data_Table const &);

private:
  const_iterator const d_begin;
  const_iterator const d_end;

  /*! We hold a copy of the scalar to prevent the problems that would arise
     *  if you took a pointer to a function-return temporary. */
  T const d_value;
};

/*! 
  Copy constructor, but update the pointers to point to the local d_value if
  they pointed to the d_value in the rhs. 
*/
template <typename T>
Data_Table<T>::Data_Table(Data_Table<T> const &rhs)
    : d_begin(rhs.d_begin), d_end(rhs.d_end), d_value(rhs.d_value) {
  if (rhs.d_begin == &(rhs.d_value)) {
    const_cast<const_iterator &>(d_begin) = &d_value;
    const_cast<const_iterator &>(d_end) = d_begin + 1;
  }
}

template <typename T>
Data_Table<T> &Data_Table<T>::operator=(Data_Table<T> const &rhs) {
  if (&rhs != this) {
    if (rhs.d_begin == &(rhs.d_value)) {
      const_cast<const_iterator &>(d_begin) = &d_value;
      const_cast<const_iterator &>(d_end) = d_begin + 1;
      const_cast<T &>(d_value) = rhs.d_value;
    } else {
      const_cast<const_iterator &>(d_begin) = rhs.d_begin;
      const_cast<const_iterator &>(d_end) = rhs.d_end;
    }
  }
  return *this;
}

/*!
  Bad things will happen if you alter the size of the source vector while
  this Data_Table is in existence.
  
  \bug Removed ctor because it does not conform to the C++ standard
  (dereferencing of end iterator is not allowed).  This particular ctor causes
  run time failures for STLport and MSVC/Debug.
*/
// template<typename T> inline
// Data_Table<T>::Data_Table(std::vector<T> const & v)
//     : d_begin(&(*(v.begin())))
//     , d_end(&(*(v.end())))
//     , d_value()
// {
//     Require(size() == v.size());
// }

template <typename T>
inline Data_Table<T>::Data_Table(const_iterator const begin,
                                 const_iterator const end)
    : d_begin(begin), d_end(end), d_value() {
  Require(!(begin > end));
}

/*! 
  Copy the scalar into a local variable, and set the pointers to that copy.
*/
template <typename T>
inline Data_Table<T>::Data_Table(T const &value)
    : d_begin(&d_value), d_end(&d_value + 1), d_value(value) {}

template <typename T>
inline Data_Table<T>::Data_Table() : d_begin(0), d_end(0), d_value() {}

template <typename T>
inline T const &Data_Table<T>::operator[](unsigned const i) const {
  Require(static_cast<int>(i) < (d_end - d_begin));
  return d_begin[i];
}

template <typename T> inline T const &Data_Table<T>::front() const {
  Require((d_end - d_begin) > 0);
  return *d_begin;
}

template <typename T> inline T const &Data_Table<T>::back() const {
  Require((d_end - d_begin) > 0);
  return *(d_end - 1);
}

template <typename T> inline T *Data_Table<T>::access() {
  Require((d_end - d_begin) > 0);
  return const_cast<T *>(d_begin);
}
}

#endif

//---------------------------------------------------------------------------//
// end of Data_Table.hh
//---------------------------------------------------------------------------//
