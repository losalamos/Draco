//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   roots/cubic1.i.hh
 * \author Kent Budge
 * \date   Wed Sep 15 10:04:02 MDT 2010
 * \brief  Solve a cubic equation assumed to have one real root
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#ifndef roots_cubic1_i_hh
#define roots_cubic1_i_hh

#include "cubic1.hh"
#include "ds++/Assert.hh"
#include "ds++/Field_Traits.hh"

namespace rtt_roots {

//---------------------------------------------------------------------------//
template <class Field>
Field cubic1(Field const &a, Field const &b, Field const &c) {
  using rtt_dsxx::value;

  double const THIRD = 1.0 / 3.0;

  Field const R = (1. / 27.) * a * a * a - (1. / 6.) * a * b + 0.5 * c;
  Field const Q = (THIRD * THIRD) * (a * a - 3.0 * b);
  Field det = R * R - Q * Q * Q;

  Check(value(R) < 0.0);
  Check(value(det) >= 0.0);
  // Assumption of single real root violated

  Field const A = pow(-R + sqrt(det), THIRD);

  Field Result = A + Q / A - a * THIRD;
  return Result;
}

} // end namespace rtt_roots

#endif // roots_cubic1_i_hh

//---------------------------------------------------------------------------//
// end of cubic1.i.hh
//---------------------------------------------------------------------------//
