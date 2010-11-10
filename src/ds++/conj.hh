//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/conj.hh
 * \author Kent Budge
 * \date   Wed Aug 11 16:04:49 2004
 * \brief  Conjugate of various field types
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_conj_hh
#define dsxx_conj_hh

#include <complex>

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the conjugate of a quantity.
 *
 * The default implementation assumes a field type that is self-conjugate,
 * such as \c double.  An example of a field type that is \em not
 * self-conjugate is \c complex.
 *
 * \arg \a Field type
 */

template<class Field>
inline Field conj(const Field &x)
{
    return x;
}

/* Specializations for non-self-conjugate types */

template<>
inline std::complex<double> conj(const std::complex<double> &x)
{
    return std::conj(x);
}

} // end namespace rtt_dsxx

#endif // dsxx_conj_hh

//---------------------------------------------------------------------------//
//              end of ds++/conj.hh
//---------------------------------------------------------------------------//
