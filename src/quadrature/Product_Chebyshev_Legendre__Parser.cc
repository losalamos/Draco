//----------------------------------*-C++-*----------------------------------------------//
/*!
 * \file   quadrature/Product_Chebyshev_Legendre.cc
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  A class representing an product Gauss-Legendre quadrature set.
 * \note   Copyright 2016-2017 Los Alamos National Security, LLC. All rights
 *         reserved. 
 */
//---------------------------------------------------------------------------------------//
// $Id: Quadrature.hh 6718 2012-08-30 20:03:01Z warsa $
//---------------------------------------------------------------------------------------//

#include <iostream>

#include "Product_Chebyshev_Legendre.hh"
#include "parser/utilities.hh"

namespace rtt_quadrature {
using namespace rtt_parser;

//---------------------------------------------------------------------------------------//
/*static*/
SP<Quadrature> Product_Chebyshev_Legendre::parse(Token_Stream &tokens) {
  // Takes two numbers, first the number of Gauss-Legendre points, which is the SN order
  Token token = tokens.shift();
  tokens.check_syntax(token.text() == "order", "expected an order");

  unsigned sn_order = parse_positive_integer(tokens);
  tokens.check_semantics(sn_order % 2 == 0, "order must be even");

  // The second number is the number of azimuthal points on each level
  // corresponding to the Gauss-Legendre points

  //std::cout << " found sn order = " << sn_order << std::endl;

  unsigned azimuthal_order = parse_positive_integer(tokens);
  tokens.check_semantics(azimuthal_order > 0,
                         "order must be greater than zero");

  //std::cout << " found azimuthal order = " << azimuthal_order << std::endl;

  bool has_axis_assignments;
  unsigned mu_axis, eta_axis;
  Octant_Quadrature::parse(tokens, has_axis_assignments, mu_axis, eta_axis);

  if (has_axis_assignments)
    return SP<Quadrature>(new Product_Chebyshev_Legendre(
        sn_order, azimuthal_order, mu_axis, eta_axis));
  else
    return SP<Quadrature>(
        new Product_Chebyshev_Legendre(sn_order, azimuthal_order));
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------------------//
//                       end of quadrature/Quadrature.cc
//---------------------------------------------------------------------------------------//
