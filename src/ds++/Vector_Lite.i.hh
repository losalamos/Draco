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
