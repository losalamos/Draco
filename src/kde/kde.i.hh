//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.i.hh
 * \author Mathew Cleveland
 * \date   Nov. 10th 2020
 * \brief  Member definitions of class kde
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef kde_kde_i_hh
#define kde_kde_i_hh

#include "ds++/dbc.hh"

namespace rtt_kde {

//------------------------------------------------------------------------------------------------//
/*!
 * reconstruction
 * \brief
 * DEFAULT reconstruction to return error if instantiation is not found
 *
 * \tparam coord enumerated value to specify KDE coordinate system
 * \tparam dim integer specifying the data dimensionality 
 * \param[in] distribution 
 * \param[in] position
 * \param[in] one_over_band_width
 * \param[inout] result returned final local function distribution
 *
 */

template <int coord>
template <int dim>
std::vector<double>
kde<coord>::reconstruction(const std::vector<double> & /*distribution*/,
                           const std::vector<std::array<double, 3>> & /*position*/,
                           const std::vector<std::array<double, 3>> & /*one_over_band_width*/,
                           const bool /*domain_decomposed*/) const {

  Insist(false, "kde::reconstruction has not been implemented for this coordinate system and or "
                "dimension combination");
  return std::vector<double>(1, 0.0);
}

//------------------------------------------------------------------------------------------------//
/*!
 * epan_kernel basis function used during reconstruction
 * \brief
 *
 * Epanechnikov kenrel to be used in reconstrtuction
 *
 * \param[in] x from kernel origin
 * \return distribution weight based on distance from the kernel center 
 *
 * Test of kde.
 */
template <int coord> double kde<coord>::epan_kernel(const double x) const {
  const double x2 = x * x;
  return x2 > 1.0 ? 0.0 : 0.75 * (1.0 - x2);
}

} // end namespace  kde

#endif // kde_kde_i_hh

//------------------------------------------------------------------------------------------------//
// end of <pkg>/kde.i.hh
//------------------------------------------------------------------------------------------------//
