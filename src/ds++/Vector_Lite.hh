//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Vector_Lite.hh
 * \brief  Header for Vector_Lite.
 * \author Rob Lowrie
 * \date   2002-10-19
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_ds_Vector_Lite_hh
#define rtt_ds_Vector_Lite_hh

#include <iostream>
#include <numeric>
#include "Assert.hh"

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
  intializers) and disadvantages (public data, cannot be a base class).
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
    inline Vector_Lite(const T &u = T());

    // Copy constructor.
    inline Vector_Lite(const Vector_Lite &rhs);

    // Constructor for N = 2.
    inline Vector_Lite(const T &u0, const T &u1);

    // Constructor for N = 3.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2);

    // Constructor for N = 4.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2, const T &u3);

    // Constructor for N = 5.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2,
		       const T &u3, const T &u4);

    /// Destructor.
    ~Vector_Lite(void) { }

    // MANIPULATORS

    // Assignment to another Vector_Lite.
    inline Vector_Lite &operator=(const Vector_Lite &rhs);

    // Assignment to a scalar.
    inline Vector_Lite &operator=(const T &rhs);

    // Comparison to another Vector_Lite.
    inline bool operator==(const Vector_Lite &a) const;

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
    bool valid_index(const int i) const
    {
	return (i >= 0 && static_cast<size_type>(i) < N);
    }

    // ...overloaded versions.
    // We need to overload to avoid automatic
    // casts to either int or unsigned versions.
    bool valid_index(const unsigned int i) const { return (i < N); }

    bool valid_index(const long int i) const { return (i >= 0 && i < N); }

    bool valid_index(const unsigned long int i) const { return (i < N); }

    // ACCESSORS

    // Indexing for int argument

    reference operator()(const int i)
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator()(const int i) const
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    reference operator[](const int i)
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator[](const int i) const
    {
	Require(valid_index(i)); return d_U[i];
    }

    // Indexing for unsigned int argument

    reference operator()(const unsigned int i)
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator()(const unsigned int i) const
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    reference operator[](const unsigned int i)
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator[](const unsigned int i) const
    {
	Require(valid_index(i)); return d_U[i];
    }

    // Indexing for long int argument

    reference operator()(const long int i)
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator()(const long int i) const
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    reference operator[](const long int i)
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator[](const long int i) const
    {
	Require(valid_index(i)); return d_U[i];
    }

    // Indexing for unsigned long int argument

    reference operator()(const unsigned long int i)
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator()(const unsigned long int i) const
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    reference operator[](const unsigned long int i)
    {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator[](const unsigned long int i) const
    {
	Require(valid_index(i)); return d_U[i];
    }

    // Iterator support
    
    iterator begin() { return d_U; }

    const_iterator begin() const { return d_U; }

    iterator end() { return d_U + N; }

    const_iterator end() const { return d_U + N; }

    size_type size() const { return N; }

    size_type max_size() const { return N; }

    bool empty() const { return N == 0; }
};

} // namespace rtt_dsxx

#include "Vector_Lite.i.hh"

#endif // rtt_ds_Vector_Lite_hh

//---------------------------------------------------------------------------//
// end of Vector_Lite.hh
//---------------------------------------------------------------------------//
