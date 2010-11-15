//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   linear/btridag.hh
 * \author Kent Budge
 * \date   Wed Sep 15 13:03:41 MDT 2010
 * \brief  Implementation of block tridiagonal solver
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef linear_btridag_hh
#define linear_btridag_hh

namespace rtt_linear
{

//---------------------------------------------------------------------------//
/*!
 * Solve a block tridiagonal system of linear equations.
 *
 * \arg \a FieldVector A random access container whose element type must
 * represent a field, such as the reals or the complex numbers. 
 *
 * \param a Subdiagonal of coefficient matrix.
 * \param b Diagonal of coefficient matrix.
 * \param c Superdiagonal of coefficient matrix.
 * \param r Right-hand side of system of equations.
 * \param n Number of blocks
 * \param m Number of variables per block
 * \param u On return, contains the solution of the system.
 *
 * \throw std::range_error If the system is not diagonal dominant.
 */
template<class FieldVector>
void btridag(FieldVector const &a,
	     FieldVector const &b,
	     FieldVector const &c,
	     FieldVector const &r,
	     unsigned const n,
	     unsigned const m,
	     FieldVector &u);

} // end namespace rtt_linear

#endif // linear_btridag_i_hh

//---------------------------------------------------------------------------//
//                 end of btridag.i.hh
//---------------------------------------------------------------------------//
