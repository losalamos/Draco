//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/ludcmp.hh
 * \author Kent Budge
 * \date   Thu Jul  1 10:54:20 2004
 * \brief  LU decomposition
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef linear_ludcmp_hh
#define linear_ludcmp_hh

namespace rtt_linear {

//---------------------------------------------------------------------------//
/*!
 * \brief LU-decompose a nonsingular matrix.
 *
 * \arg \a FieldVector1 A random-access container type on a field.
 * \arg \a IntVector A random-access container type on an integral type.
 *
 * \param a Matrix to decompose.  On return, contains the decomposition.
 * \param indx On return, contains the pivoting map.
 * \param d On return, contains the sign of the determinant.
 *
 * \pre \c a.size()==indx.size()*indx.size()
 */
template <class FieldVector, class IntVector>
void ludcmp(FieldVector &a, IntVector &indx,
            typename FieldVector::value_type &d);

//---------------------------------------------------------------------------//
/*!
 * \brief Solve the system \f$Ax=b\f$
 *
 * \arg \a FieldVector1 A random-access container type on a field.
 * \arg \a IntVector A random-access container type on an integral type.
 * \arg \a FieldVector2 A random-access container type on a field.
 *
 * \param a LU decomposition of \f$A\f$.
 * \param indx Pivot map for decomposition of \f$A\f$.
 * \param b Right-hand side \f$b\f$.  On return, contains solution \f$x\f$.
 *
 * \pre \c a.size()==indx.size()*indx.size()
 * \pre \c b.size()==indx.size()
 */
template <class FieldVector1, class IntVector, class FieldVector2>
void lubksb(FieldVector1 const &a, IntVector const &indx, FieldVector2 &b);

} // end namespace astr_star1d

#endif // linear_ludcmp_hh

//---------------------------------------------------------------------------//
// end of Implicit.cc
//---------------------------------------------------------------------------//
