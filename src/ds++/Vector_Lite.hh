//----------------------------------*-C++-*----------------------------------//
/*!
  \file   Vector_Lite.hh
  \brief  Header for Vector_Lite.
  \author Rob Lowrie
  \date   2002-10-19
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

  It adds iterator and arithemtic support, along with bounds checking
  (via Draco's DBC).

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

    /// Constructor based on a single element type.
    /// Sets all values to \a u.
    Vector_Lite(const T &u = T()) {
	for ( size_type i = 0; i < N; i++ ) {
	    d_U[i] = u;
	}
    }

    /// Copy constructor.
    Vector_Lite(const Vector_Lite &rhs) {
	for ( size_type i = 0; i < N; i++ ) {
	    d_U[i] = rhs.d_U[i];
	}
    }

    /// Constructor for N = 2.
    Vector_Lite(const T &u0, const T &u1) {
 	Require(N == 2); d_U[0] = u0; d_U[1] = u1;
    }

    /// Constructor for N = 3.
    Vector_Lite(const T &u0, const T &u1, const T &u2) {
	Require(N == 3); d_U[0] = u0; d_U[1] = u1; d_U[2] = u2;
    }

    /// Constructor for N = 4.
    Vector_Lite(const T &u0, const T &u1, const T &u2,
		const T &u3) {
	Require(N == 4); d_U[0] = u0; d_U[1] = u1; d_U[2] = u2;
	d_U[3] = u3;
    }

    /// Constructor for N = 5.
    Vector_Lite(const T &u0, const T &u1, const T &u2,
	    const T &u3, const T &u4) {
	Require(N == 5); d_U[0] = u0; d_U[1] = u1; d_U[2] = u2;
	d_U[3] = u3; d_U[4] = u4;
    }

    /// Constructor for N = 6.
    Vector_Lite(const T &u0, const T &u1, const T &u2,
		const T &u3, const T &u4, const T &u5) {
	Require(N == 6); d_U[0] = u0; d_U[1] = u1; d_U[2] = u2;
	d_U[3] = u3; d_U[4] = u4; d_U[5] = u5;
    }

    /// Destructor.
    virtual ~Vector_Lite(void) { }

    // MANIPULATORS

    /// Assignment to another Vector_Lite.
    Vector_Lite &operator=(const Vector_Lite &rhs) {
	if ( this == &rhs ) {
	    return *this;
	}
	
	for ( size_type i = 0; i < N; i++ )
	    d_U[i] = rhs.d_U[i];
	
	return *this;
    }

    /// Assignment to type T.
    Vector_Lite &operator=(const T rhs) {
	
	for ( size_type i = 0; i < N; i++ )
	    d_U[i] = rhs;

	return *this;
    }

    /// Comparison to another Vector_Lite.
    bool operator==(const Vector_Lite &a) const {
	for ( size_type i = 0; i < N; i++ )
	    if ( d_U[i] != a.d_U[i] ) return false;
	return true;
    }

    // Basic arithmetic operations

    Vector_Lite &operator+=(const Vector_Lite &a) {
	for ( size_type i = 0; i < N; i++ )
	    d_U[i] += a.d_U[i];
	return *this;
    }

    Vector_Lite &operator-=(const Vector_Lite &a) {
	for ( size_type i = 0; i < N; i++ )
	    d_U[i] -= a.d_U[i];
	return *this;
    }

    Vector_Lite &operator*=(const Vector_Lite &a) {
	for ( size_type i = 0; i < N; i++ )
	    d_U[i] *= a.d_U[i];
	return *this;
    }

    Vector_Lite &operator/=(const Vector_Lite &a) {
	for ( size_type i = 0; i < N; i++ )
	    d_U[i] /= a.d_U[i];
	return *this;
    }

    Vector_Lite &operator+=(const T a) {
	for ( size_type i = 0; i < N; i++ )
	    d_U[i] += a;
	return *this;
    }

    Vector_Lite &operator-=(const T a) {
	for ( size_type i = 0; i < N; i++ )
	    d_U[i] -= a;
	return *this;
    }

    Vector_Lite &operator*=(const T a) {
	for ( size_type i = 0; i < N; i++ )
	    d_U[i] *= a;
	return *this;
    }

    Vector_Lite &operator/=(const T a) {
	for ( size_type i = 0; i < N; i++ )
	    d_U[i] /= a;
	return *this;
    }

    /// Returns true if \a i is a valid array index.
    bool valid_index(const int i) const {
	return (i >= 0 && static_cast<size_type>(i) < N);
    }

    // ...overloaded versions.
    // We need to overload to avoid automatic
    // casts to either int or unsigned versions.
    bool valid_index(const unsigned int i) const {
	return (i < N);
    }

    bool valid_index(const long int i) const {
	return (i >= 0 && i < N);
    }

    bool valid_index(const unsigned long int i) const {
	return (i < N);
    }

    // ACCESSORS

    // Indexing for int argument

    reference operator()(const int i) {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator()(const int i) const {
	Require(valid_index(i)); return d_U[i];
    }
    
    reference operator[](const int i) {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator[](const int i) const {
	Require(valid_index(i)); return d_U[i];
    }

    // Indexing for unsigned int argument

    reference operator()(const unsigned int i) {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator()(const unsigned int i) const {
	Require(valid_index(i)); return d_U[i];
    }
    
    reference operator[](const unsigned int i) {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator[](const unsigned int i) const {
	Require(valid_index(i)); return d_U[i];
    }

    // Indexing for long int argument

    reference operator()(const long int i) {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator()(const long int i) const {
	Require(valid_index(i)); return d_U[i];
    }
    
    reference operator[](const long int i) {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator[](const long int i) const {
	Require(valid_index(i)); return d_U[i];
    }

    // Indexing for unsigned long int argument

    reference operator()(const unsigned long int i) {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator()(const unsigned long int i) const {
	Require(valid_index(i)); return d_U[i];
    }
    
    reference operator[](const unsigned long int i) {
	Require(valid_index(i)); return d_U[i];
    }
    
    const_reference operator[](const unsigned long int i) const {
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

//---------------------------------------------------------------------------//
// free functions
//---------------------------------------------------------------------------//

template <class T, size_t N>
T inner_product(const rtt_dsxx::Vector_Lite<T, N> &a,
		const rtt_dsxx::Vector_Lite<T, N> &b)
{
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

} // namespace rtt_dsxx

//---------------------------------------------------------------------------//
// global opeators
//---------------------------------------------------------------------------//

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator+(const rtt_dsxx::Vector_Lite<T, N> &a,
	  const rtt_dsxx::Vector_Lite<T, N> &b)
{
    return rtt_dsxx::Vector_Lite<T, N>(a) += b;
}

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator-(const rtt_dsxx::Vector_Lite<T, N> &a,
	  const rtt_dsxx::Vector_Lite<T, N> &b)
{
    return rtt_dsxx::Vector_Lite<T, N>(a) -= b;
}

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator*(const rtt_dsxx::Vector_Lite<T, N> &a,
	  const rtt_dsxx::Vector_Lite<T, N> &b)
{
    return rtt_dsxx::Vector_Lite<T, N>(a) *= b;
}

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator/(const rtt_dsxx::Vector_Lite<T, N> &a,
	  const rtt_dsxx::Vector_Lite<T, N> &b)
{
    return rtt_dsxx::Vector_Lite<T, N>(a) /= b;
}

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator+(const rtt_dsxx::Vector_Lite<T, N> &a,
	  const T b)
{
    return rtt_dsxx::Vector_Lite<T, N>(a) += b;
}

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator+(const T a,
	  const rtt_dsxx::Vector_Lite<T, N> &b)
{
    return rtt_dsxx::Vector_Lite<T, N>(b) += a;
}

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator-(const rtt_dsxx::Vector_Lite<T, N> &a,
	  const T b)
{
    return rtt_dsxx::Vector_Lite<T, N>(a) -= b;
}

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator-(const T a,
	  const rtt_dsxx::Vector_Lite<T, N> &b)
{
    return rtt_dsxx::Vector_Lite<T, N>(b) -= a;
}

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator*(const rtt_dsxx::Vector_Lite<T, N> &a,
	  const T b)
{
    return rtt_dsxx::Vector_Lite<T, N>(a) *= b;
}

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator*(const T a,
	  const rtt_dsxx::Vector_Lite<T, N> &b)
{
    return rtt_dsxx::Vector_Lite<T, N>(b) *= a;
}

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator/(const rtt_dsxx::Vector_Lite<T, N> &a,
	  const T b)
{
    return rtt_dsxx::Vector_Lite<T, N>(a) /= b;
}

template <class T, size_t N>
inline const rtt_dsxx::Vector_Lite<T, N>
operator-(const rtt_dsxx::Vector_Lite<T, N> &a)
{
    rtt_dsxx::Vector_Lite<T, N> neg(a);
    
    for (size_t i = 0; i < N; i++ )
	neg(i) = -a(i);
    
    return neg;
}

//---------------------------------------------------------------------------//
// stream opeators
//---------------------------------------------------------------------------//

template <class T, size_t N>
std::ostream &operator<<(std::ostream &os,
			 const rtt_dsxx::Vector_Lite<T, N> &a)
{
    os << a(0);
    
    for (size_t i = 1; i < N; i++) {
	os << " " << a(i);
    }
    
    return os;
}

template <class T, size_t N>
std::istream &operator>>(std::istream &os,
			 rtt_dsxx::Vector_Lite<T, N> &a)
{
    for (size_t i = 0; i < N; i++) {
	os >> a(i);
    }
    
    return os;
}

#endif // rtt_ds_Vector_Lite_hh

//---------------------------------------------------------------------------//
// end of Vector_Lite.hh
//---------------------------------------------------------------------------//
