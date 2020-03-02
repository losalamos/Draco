//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   quadrature/Lobatto.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an interval Gauss-Legendre quadrature set.
 * \note   Copyright 2016-2020 Triad National Security, LLC. All rights
 *         reserved. */
//----------------------------------------------------------------------------//

#ifndef rtt_quadrature_Lobatto_hh
#define rtt_quadrature_Lobatto_hh

#include "Interval_Quadrature.hh"

namespace rtt_quadrature {

//============================================================================//
/*!
 * \class Lobatto
 *
 * \brief A class representing an interval Lobatto quadrature set.
 *
 * This is an interval (e.g. 1D) angle quadrature set whose abscissae at order N
 * are the roots of the derivative of the Laguerre polynomial of order N-1 plus
 * the end points of the interval.
 */
//============================================================================//

class Lobatto : public Interval_Quadrature {
public:
  // CREATORS
  DLL_PUBLIC_quadrature explicit Lobatto(unsigned sn_order);

  // ACCESSORS

  // SERVICES

  virtual string name() const;

  virtual string parse_name() const;

  virtual unsigned number_of_levels() const;

  virtual string as_text(string const &indent) const;

  virtual bool is_open_interval() const;

  bool check_class_invariants() const;

  // STATICS

  static std::shared_ptr<Quadrature> parse(Token_Stream &tokens);

protected:
  virtual vector<Ordinate> create_level_ordinates_(double norm) const;
};

} // end namespace rtt_quadrature

#endif // rtt_quadrature_Quadrature_hh

//----------------------------------------------------------------------------//
// end of quadrature/Quadrature.hh
//----------------------------------------------------------------------------//
