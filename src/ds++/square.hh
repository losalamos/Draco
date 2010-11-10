//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/square.hh
 * \author Kent Budge
 * \date   Tue Jul  6 08:57:41 2004
 * \brief  Declare the square function template.
 * \note   © Copyright 2006 LANSLLC All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef dsxx_square_hh
#define dsxx_square_hh

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the square of a value.
 * 
 * \arg \a Semigroup A type representing an algebraic structure closed under
 * multiplication, such as the integers or the reals.
 *
 * \param x Value to be squared.
 * \return \f$x^2\f$
 */
template <class Semigroup>
inline Semigroup square(const Semigroup &x)
{
    return x*x;
}

} // end namespace rtt_dsxx

#endif // dsxx_square_hh

//---------------------------------------------------------------------------//
//              end of dsxx/square.hh
//---------------------------------------------------------------------------//
