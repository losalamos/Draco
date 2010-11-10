//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/cube.hh
 * \author Kent Budge
 * \date   Tue Jul  6 10:03:25 MDT 2004
 * \brief  Declare the cube function template.
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_cube_hh
#define dsxx_cube_hh

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the cube of a value.
 * 
 * \arg \a Semigroup A type representing an algebraic structure closed under
 * multiplication such as the integers or the reals.
 *
 * \param x Value to be cubed.
 * \return \f$x^3\f$
 */
template <class Semigroup>
inline Semigroup cube(Semigroup const &x)
{
    return x*x*x;
}

} // end namespace rtt_dsxx

#endif // dsxx_cube_hh

//---------------------------------------------------------------------------//
//              end of ds++/cube.hh
//---------------------------------------------------------------------------//
