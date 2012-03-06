//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Vector_Lite.hh
 * \brief  Header for Vector_Lite.
 * \note   Copyright (C) 2003-2012 Los Alamos National Security, LLC
 * \version $Id$
 */
//---------------------------------------------------------------------------//

#ifndef rtt_ds_Vector_Lite_hh
#define rtt_ds_Vector_Lite_hh

#include "Assert.hh"
#include <cstddef>
#include <iostream>
#include <numeric>

namespace rtt_dsxx
{

//===========================================================================//
/*!
  \class Vector_Lite

  \brief  Array container that is a wrapper around a standard C array.

  It adds iterator and arithemtic support, along with bounds checking (via
  Draco's DBC).

  An alternative to this class is boost::array (www.boost.org).  However,
  boost::array is an aggregate type, which has advantages (can use
  initializers) and disadvantages (public data, cannot be a base class).
  boost::array also doesn't do bounds checking.
  
  \param T Type of each array element.
  
  \param N Length of array.
*/
/*!
 * \example ds++/test/tstVector_Lite.cc
 *
 * Test of Vector_Lite.
 */
//===========================================================================//
template <class T, size_t N>
class Vector_Lite
{
  public:

    // TYPEDEFS AND DATA TYPES

    typedef T value_type;
    typedef value_type* pointer;
    typedef const T* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef ptrdiff_t difference_type;
    typedef size_t size_type;
    typedef pointer iterator;
    typedef const_pointer const_iterator;
    
  private:

    // DATA

    T d_U[N];
    
  public:

    // CREATORS

    // Constructor based on a scalar value.
    inline explicit Vector_Lite(const T &u = T());

    // Constructor for N = 2.
    inline Vector_Lite(const T &u0, const T &u1);

    // Constructor for N = 3.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2);

    // Constructor for N = 4.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2, const T &u3);

    // Constructor for N = 5.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2,
		       const T &u3, const T &u4);

    // Fill in from C array.
    inline void fill(const T u[N]);

    /// Destructor.
    ~Vector_Lite(void) { }

    // MANIPULATORS

    // Assignment to a another Vector_Lite.
    inline Vector_Lite &operator=(const Vector_Lite &rhs);

    // Assignment to a scalar.
    inline Vector_Lite &operator=(const T &rhs);

    // Comparisons to another Vector_Lite.
    inline bool operator==(const Vector_Lite &a) const;
    inline bool operator<(const Vector_Lite &a) const;

    // Basic arithmetic operations, vector right-hand side

    inline Vector_Lite &operator+=(const Vector_Lite &a);
    inline Vector_Lite &operator-=(const Vector_Lite &a);
    inline Vector_Lite &operator*=(const Vector_Lite &a);
    inline Vector_Lite &operator/=(const Vector_Lite &a);

    // Basic arithmetic operations, scalar right-hand side

    inline Vector_Lite &operator+=(const T &a);
    inline Vector_Lite &operator-=(const T &a);
    inline Vector_Lite &operator*=(const T &a);
    inline Vector_Lite &operator/=(const T &a);
    
    /// Returns true if \a i is a valid array index.
    bool valid_index(const size_type i) const
    {
	return i < N;
    }


    // ACCESSORS

    // Indexing for int argument

    /// Indexing using ()
    reference operator()(const size_type i)
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator()(const size_type i) const
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    /// Indexing using []
    reference operator[](const size_type i)
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    /// const indexing using []
    const_reference operator[](const size_type i) const
    {
	Require(valid_index(i)); return d_U[i];
    }

    // Iterator support
    
    /// iterator begin
    iterator begin() { return d_U; }

    /// const iterator begin
    const_iterator begin() const { return d_U; }

    /// iterator end
    iterator end() { return d_U + N; }

    /// const iterator end
    const_iterator end() const { return d_U + N; }

    /// number of elements (\a N); for STL support
    size_type size() const { return N; }

    /// max number of elements (\a N); for STL support
    size_type max_size() const { return N; }

    /// true if \a N = 0; for STL support
    bool empty() const { return N == 0; }
};

} // namespace rtt_dsxx

#include "Vector_Lite.i.hh"

#endif // rtt_ds_Vector_Lite_hh

//---------------------------------------------------------------------------//
// end of Vector_Lite.hh
//---------------------------------------------------------------------------//
