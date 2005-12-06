//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/test/cube.hh
 * \author Kent Budge
 * \date   Tue Jul  6 10:03:25 MDT 2004
 * \brief  Declare the cube function template.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_cube_hh
#define quadrature_cube_hh

namespace rtt_quadrature
{

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the cube of a value.
 * 
 * \arg \a Group A type representing a mathematical group, such as the
 * integers or the reals.
 * \param x Value to be cubed.
 * \return \f$x^3\f$
 */
template <class Group>
inline Group cube(Group const &x)
{
    return x*x*x;
}

} // end namespace rtt_quadrature

#endif // quadrature_cube_hh

//---------------------------------------------------------------------------//
//              end of quadrature/test/cube.hh
//---------------------------------------------------------------------------//
