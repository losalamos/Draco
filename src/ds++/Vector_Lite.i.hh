//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/Vector_Lite.i.hh
 * \brief  Implementation of functions for Vector_Lite.
 * \author Rob Lowrie
 * \date   2002-10-19
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_ds_Vector_Lite_i_hh
#define rtt_ds_Vector_Lite_i_hh

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
// member functions
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor based on a scalar value.
 *
 * Initializes all values to \a u.  Note that this ctor also acts as the
 * default ctor.
 * 
 * \param u     Scalar value.
 */
template <class T, size_t N>
Vector_Lite<T, N>::
Vector_Lite(const T &u)
{
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] = u;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Copy constructor.
 * 
 * \param rhs   Value to copy.
 */
template <class T, size_t N>
Vector_Lite<T, N>::
Vector_Lite(const Vector_Lite &rhs)
{
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] = rhs.d_U[i];
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 2.
 * 
 * \param u0 1st element.    
 * \param u1 2nd element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::
Vector_Lite(const T &u0,
	    const T &u1)
{
    Require(N == 2);
    d_U[0] = u0;
    d_U[1] = u1;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 3.
 * 
 * \param u0 1st element.    
 * \param u1 2nd element.
 * \param u2 3rd element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::
Vector_Lite(const T &u0,
	    const T &u1,
	    const T &u2)
{
    Require(N == 3);
    d_U[0] = u0;
    d_U[1] = u1;
    d_U[2] = u2;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 4.
 * 
 * \param u0 1st element.    
 * \param u1 2nd element.
 * \param u2 3rd element.
 * \param u3 4th element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::
Vector_Lite(const T &u0,
	    const T &u1,
	    const T &u2,
	    const T &u3)
{
    Require(N == 4);
    d_U[0] = u0;
    d_U[1] = u1;
    d_U[2] = u2;
    d_U[3] = u3;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 5.
 * 
 * \param u0 1st element.    
 * \param u1 2nd element.
 * \param u2 3rd element.
 * \param u3 4th element.
 * \param u4 5th element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::
Vector_Lite(const T &u0,
	    const T &u1,
	    const T &u2,
	    const T &u3,
	    const T &u4)
{
    Require(N == 5);
    d_U[0] = u0;
    d_U[1] = u1;
    d_U[2] = u2;
    d_U[3] = u3;
    d_U[4] = u4;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment to another Vector_Lite.
 */
template <class T, size_t N>
Vector_Lite<T, N> &
Vector_Lite<T, N>::
operator=(const Vector_Lite<T, N> &rhs)
{
    if ( this == &rhs )
    {
	return *this;
    }
	
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] = rhs.d_U[i];
    }
	
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment to a scalar.
 */
template <class T, size_t N>
Vector_Lite<T, N> &
Vector_Lite<T, N>::
operator=(const T &rhs)
{
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] = rhs;
    }
	
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Comparison to another Vector_Lite.
 */
template <class T, size_t N>
bool
Vector_Lite<T, N>::
operator==(const Vector_Lite<T, N> &a) const
{
    for ( size_type i = 0; i < N; i++ )
    {
	if ( d_U[i] != a.d_U[i] ) return false;
    }
    
    return true;
}

//---------------------------------------------------------------------------//
// basic arithmetic member functions, vector right-hand side
//---------------------------------------------------------------------------//

template <class T, size_t N>
Vector_Lite<T,N> &
Vector_Lite<T, N>::
operator+=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] += a.d_U[i];
    }
    
    return *this;
}

template <class T, size_t N>
Vector_Lite<T,N> &
Vector_Lite<T, N>::
operator-=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] -= a.d_U[i];
    }
    
    return *this;
}

template <class T, size_t N>
Vector_Lite<T,N> &
Vector_Lite<T, N>::
operator*=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] *= a.d_U[i];
    }
    
    return *this;
}

template <class T, size_t N>
Vector_Lite<T,N> &
Vector_Lite<T, N>::
operator/=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] /= a.d_U[i];
    }
    
    return *this;
}

//---------------------------------------------------------------------------//
// basic arithmetic member functions, scalar right-hand side
//---------------------------------------------------------------------------//

template <class T, size_t N>
Vector_Lite<T,N> &
Vector_Lite<T, N>::
operator+=(const T &a)
{
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] += a;
    }
    
    return *this;
}

template <class T, size_t N>
Vector_Lite<T,N> &
Vector_Lite<T, N>::
operator-=(const T &a)
{
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] -= a;
    }
    
    return *this;
}

template <class T, size_t N>
Vector_Lite<T,N> &
Vector_Lite<T, N>::
operator*=(const T &a)
{
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] *= a;
    }
    
    return *this;
}

template <class T, size_t N>
Vector_Lite<T,N> &
Vector_Lite<T, N>::
operator/=(const T &a)
{
    for ( size_type i = 0; i < N; i++ )
    {
	d_U[i] /= a;
    }
    
    return *this;
}

//---------------------------------------------------------------------------//
// free functions
//---------------------------------------------------------------------------//

template <class T, size_t N>
T inner_product(const Vector_Lite<T, N> &a,
		const Vector_Lite<T, N> &b)
{
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

//---------------------------------------------------------------------------//
// global operators
//---------------------------------------------------------------------------//

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator+(const Vector_Lite<T, N> &a,
	  const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) += b;
}

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator-(const Vector_Lite<T, N> &a,
	  const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) -= b;
}

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator*(const Vector_Lite<T, N> &a,
	  const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) *= b;
}

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator/(const Vector_Lite<T, N> &a,
	  const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) /= b;
}

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator+(const Vector_Lite<T, N> &a,
	  const T b)
{
    return Vector_Lite<T, N>(a) += b;
}

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator+(const T a,
	  const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) += a;
}

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator-(const Vector_Lite<T, N> &a,
	  const T b)
{
    return Vector_Lite<T, N>(a) -= b;
}

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator-(const T a,
	  const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) -= a;
}

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator*(const Vector_Lite<T, N> &a,
	  const T b)
{
    return Vector_Lite<T, N>(a) *= b;
}

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator*(const T a,
	  const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) *= a;
}

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator/(const Vector_Lite<T, N> &a,
	  const T b)
{
    return Vector_Lite<T, N>(a) /= b;
}

template <class T, size_t N>
inline const Vector_Lite<T, N>
operator-(const Vector_Lite<T, N> &a)
{
    Vector_Lite<T, N> neg(a);
    
    for (size_t i = 0; i < N; i++ )
	neg(i) = -a(i);
    
    return neg;
}

//---------------------------------------------------------------------------//
// stream opeators
//---------------------------------------------------------------------------//

template <class T, size_t N>
std::ostream &operator<<(std::ostream &os,
			 const Vector_Lite<T, N> &a)
{
    os << a(0);
    
    for (size_t i = 1; i < N; i++) {
	os << " " << a(i);
    }
    
    return os;
}

template <class T, size_t N>
std::istream &operator>>(std::istream &os,
			 Vector_Lite<T, N> &a)
{
    for (size_t i = 0; i < N; i++) {
	os >> a(i);
    }
    
    return os;
}

} // namespace rtt_dsxx

#endif // rtt_ds_Vector_Lite_i_hh

//---------------------------------------------------------------------------//
// end of Vector_Lite.i.hh
//---------------------------------------------------------------------------//
