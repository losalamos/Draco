//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/DracoMath.hh
 * \author Kent G. Budge
 * \date   Wed Jan 22 15:18:23 MST 2003
 * \brief  New or overloaded cmath or cmath-like functions.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_DracoMath_hh
#define rtt_dsxx_DracoMath_hh

#include "Constexpr_Functions.hh"
#include "Soft_Equivalence.hh"
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <functional>

namespace rtt_dsxx {

//---------------------------------------------------------------------------//
// isFinite.hh
//
// Try to use the C++11/C99 functions isinf, isnan and isfinite defined in
// <cmath> instead of defining our own.  I would like to use C++11
// implemenations which are true functions in the std:: namespace.  The problem
// here is that PGI/13.7 does not have these language features.  However, PGI
// does provide the C99 _macros_ of the same name (w/o namespace qualifier).
// ---------------------------------------------------------------------------//
#if defined _WIN32 || defined __CYGWIN__

template <typename T> bool isNan(T a) { return _isnan(a); }
template <typename T> bool isInf(T a) { return !_finite(a); }
template <typename T> bool isFinite(T a) { return _finite(a); }

#elif defined draco_isPGI

template <typename T> bool isNan(T a) { return isnan(a); }
template <typename T> bool isInf(T a) { return isinf(a); }
template <typename T> bool isFinite(T a) { return isfinite(a); }

#else

template <typename T> bool isNan(T a) { return std::isnan(a); }
template <typename T> bool isInf(T a) { return std::isinf(a); }
template <typename T> bool isFinite(T a) { return std::isfinite(a); }

#endif

//---------------------------------------------------------------------------//
/*!
 * \brief Return the conjugate of a quantity.
 *
 * The default implementation assumes a field type that is self-conjugate, such
 * as \c double.  An example of a field type that is \em not self-conjugate is
 * \c complex.
 *
 * \tparam Field type
 * \param arg Field type
 */
template <typename Field> inline Field conj(const Field &arg) { return arg; }

// Specializations for non-self-conjugate types
template <> inline std::complex<double> conj(const std::complex<double> &arg) {
  return std::conj(arg);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the cube of a value.
 *
 * \param[in] x Value to be cubed.
 * \return \f$x^3\f$
 *
 * \c Semigroup is a type representing an algebraic structure closed under
 * multiplication such as the integers or the reals.
 */
template <typename Semigroup> inline Semigroup cube(Semigroup const &x) {
  return x * x * x;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Return the positive difference of the arguments.
 *
 * This is a replacement for the FORTRAN DIM function.
 *
 * \arg \a Ordered_Group_Element A type for which operator< and unary operator-
 *      are defined and which can be constructed from a literal \c 0.
 *
 * \param a Minuend
 * \param b Subtrahend
 * \return \f$max(0, a-b)\f$
 *
 * \deprecated A FORTRAN relic that should disappear eventually.
 */
template <typename Ordered_Group_Element>
inline Ordered_Group_Element dim(Ordered_Group_Element a,
                                 Ordered_Group_Element b) {
  if (a < b)
    return Ordered_Group_Element(0);
  else
    return a - b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the square of a value.
 *
 * \arg \a Semigroup A type representing an algebraic structure closed under
 *      multiplication, such as the integers or the reals.
 *
 * \param[in] x Value to be squared.
 * \return \f$x^2\f$
 */
template <typename Semigroup> inline Semigroup square(const Semigroup &x) {
  return x * x;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the hypotenuse of a right triangle.
 *
 * This function evaluates the expression \f$\sqrt{a^2+b^2}\f$ in a way that is
 * insensitive to roundoff and preserves range.
 *
 * \arg \a Real A real number type
 * \param a First leg of triangle
 * \param b Second leg of triangle
 * \return Hypotenuse, \f$\sqrt{a^2+b^2}\f$
 */
template <typename Real> inline double pythag(Real a, Real b) {
  using std::abs;
  Real absa = abs(a), absb = abs(b);
  // We must avoid (a/b)^2 > max.
  if (absa <= absb * std::sqrt(std::numeric_limits<Real>::min()))
    return absb;
  if (absb <= absa * std::sqrt(std::numeric_limits<Real>::min()))
    return absa;
  // The regular case...
  if (absa > absb)
    return absa * std::sqrt(1.0 + square(absb / absa));
  else
    return absb * std::sqrt(1.0 + square(absa / absb));
}

//---------------------------------------------------------------------------//
/*!
 * \brief  Transfer the sign of the second argument to the first argument.
 *
 * This is a replacement for the FORTRAN SIGN function.  It is useful in
 * numerical algorithms that are roundoff or overflow insensitive and should not
 * be deprecated.
 *
 * \arg \a Ordered_Group
 * A type for which \c operator< and unary \c operator- are defined and which
 * can be compared to literal \c 0.
 *
 * \param a Argument supplying magnitude of result.
 * \param b Argument supplying sign of result.
 * \return \f$|a|sgn(b)\f$
 */
template <typename Ordered_Group>
inline Ordered_Group sign(Ordered_Group a, Ordered_Group b) {
  using std::abs; // just to be clear

  if (b < 0)
    return -abs(a);
  else
    return abs(a);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do a linear interpolation between two values.
 *
 * \param[in] x1 x coordinate of first data point.
 * \param[in] y1 y coordinate of first data point.
 * \param[in] x2 x coordinate of second data point.
 * \param[in] y2 y coordinate of second data point.
 * \param[in] x  x coordinate associated with requested y value.
 * \return The y value associated with x based on linear interpolation between
 *         (x1,y1) and (x2,y2).
 *
 * Given two points (x1,y1) and (x2,y2), use linaer interpolation to find the y
 * value associated with the provided x value.
 *
 *          y2-y1
 * y = y1 + ----- * (x-x1)
 *          x2-x1
 *
 * \pre  x in (x1,x2), extrapolation is not allowed.
 * \post y in (y1,y2), extrapolation is not allowed.
 */
constexpr inline double linear_interpolate(double const x1, double const x2,
                                           double const y1, double const y2,
                                           double const x) {
  Require(ce_fabs(x2 - x1) > std::numeric_limits<double>::epsilon());
  Require(((x >= x1) && (x <= x2)) || ((x >= x2) && (x <= x1)));

  // return value
  double const value = (y2 - y1) / (x2 - x1) * (x - x1) + y1;

  Ensure(((value >= y1) && (value <= y2)) || ((value >= y2) && (value <= y1)));
  return value;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do a 3D linear interpolation between vertices of a rectangular prism.
 *
 * Algorithm from wikipedia's Trilinear Interpolation article, hat tip to E.
 * Norris for the reference.
 *
 * \param[in] x0   lower x coordinate of lattice
 * \param[in] x1   upper x coordinate of lattic
 * \param[in] y0   lower y coordinate of lattice
 * \param[in] y1   upper y coordinate of lattice
 * \param[in] z0   lower z coordinate of lattice
 * \param[in] z1   upper z coordinate of lattic
 * \param[in] f000 function at (x0,y0,z0)
 * \param[in] f100 function at (x1,y0,z0)
 * \param[in] f001 function at (x0,y0,z1)
 * \param[in] f101 function at (x1,y0,z1)
 * \param[in] f010 function at (x0,y1,z0)
 * \param[in] f110 function at (x1,y1,z0)
 * \param[in] f011 function at (x0,y1,z1)
 * \param[in] f111 function at (x1,y1,z1)
 * \param[in] x    x coordinate of interpolation point
 * \param[in] y    y coordinate of interpolation point
 * \param[in] z    z coordinate of interpolation point
 * \return The function value linearly interpolated to (x,y,z)
 */
inline double
linear_interpolate_3(double const x0, double const x1, double const y0,
                     double const y1, double const z0, double const z1,
                     double const f000, double const f100, double const f001,
                     double const f101, double const f010, double const f110,
                     double const f011, double const f111, double const x,
                     double const y, double const z) {
  Require(std::abs(x1 - x0) > std::numeric_limits<double>::epsilon());
  Require(std::abs(y1 - y0) > std::numeric_limits<double>::epsilon());
  Require(std::abs(z1 - z0) > std::numeric_limits<double>::epsilon());
  Require((x >= x0) && (x <= x1) && (y >= y0) && (y <= y1) && (z >= z0) &&
          (z <= z1));

  double xd = (x - x0) / (x1 - x0);
  double yd = (y - y0) / (y1 - y0);
  double zd = (z - z0) / (z1 - z0);

  double f00 = f000 * (1. - xd) + f100 * xd;
  double f01 = f001 * (1. - xd) + f101 * xd;
  double f10 = f010 * (1. - xd) + f110 * xd;
  double f11 = f011 * (1. - xd) + f111 * xd;

  double f0 = f00 * (1. - yd) + f10 * yd;
  double f1 = f01 * (1. - yd) + f11 * yd;

  double f = f0 * (1. - zd) + f1 * zd;

  return f;
}

} // namespace rtt_dsxx

#endif // rtt_dsxx_DracoMath_hh

//---------------------------------------------------------------------------//
// end of DracoMath.hh
//---------------------------------------------------------------------------//
