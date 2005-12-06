//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/square.hh
 * \author Kent Budge
 * \date   Tue Jul  6 08:57:41 2004
 * \brief  Declare the square function template.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_square_hh
#define quadrature_square_hh

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the square of a value.
 * 
 * \arg \a Group A type representing a mathematical group, such as the
 * integers or the reals.
 * \param x Value to be squared.
 * \return \f$x^2\f$
 */
template <class Group>
inline Group square(const Group &x)
{
    return x*x;
}

} // end namespace rtt_quadrature

#endif // quadrature_square_hh

//---------------------------------------------------------------------------//
//              end of quadrature/test/square.hh
//---------------------------------------------------------------------------//
