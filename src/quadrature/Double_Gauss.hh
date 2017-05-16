//-----------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Double_Gauss.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval double Gauss-Legendre quadrature
 *         set.
 * \note   Copyright 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef __quadrature_Double_Gauss_hh__
#define __quadrature_Double_Gauss_hh__

#include "Interval_Quadrature.hh"

namespace rtt_quadrature {

//============================================================================//
/*!
 * \class Double_Gauss
 *
 * \brief A class representing an interval double Gauss-Legendre quadrature set.
 *
 * This is an interval (e.g. 1D) angle quadrature set in which the
 * Gauss-Legendre quadrature of order n/2 is mapped separately onto the
 * intervals [-1,0] and [0,1] assuming N>2. The case N=2 does not integrate the
 * flux and we substitute ordinary Gauss-Legendre quadrature (N=2) on the full
 * interval [-1,1].
 */
//============================================================================//

class Double_Gauss : public Interval_Quadrature {
public:
  // CREATORS
  DLL_PUBLIC_quadrature explicit Double_Gauss(unsigned sn_order);

  // ACCESSORS

  // SERVICES

  virtual string name() const;

  virtual string parse_name() const;

  virtual unsigned number_of_levels() const;

  virtual string as_text(string const &indent) const;

  bool check_class_invariants() const;

  // STATICS

  static std::shared_ptr<Quadrature> parse(Token_Stream &tokens);

protected:
  virtual vector<Ordinate> create_level_ordinates_(double norm) const;

  // DATA
};

} // end namespace rtt_quadrature

#endif // __quadrature_Quadrature_hh__

//----------------------------------------------------------------------------//
// end of quadrature/Quadrature.hh
//----------------------------------------------------------------------------//
