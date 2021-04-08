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
 * \brief
 * DEFAULT reconstruction to return error if instantiation is not found
 *
 * \tparam coord enumerated value to specify KDE coordinate system
 * \tparam dim integer specifying the data dimensionality 
 * \param[in] distribution 
 * \param[in] position
 * \param[in] one_over_band_width
 * \param[in] domain_decomposed
 * \return final local function distribution
 *
 */

template <int coord>
template <int dim>
std::vector<double>
kde<coord>::reconstruction(const std::vector<double> & /*distribution*/,
                           const std::vector<std::array<double, 3>> & /*position*/,
                           const std::vector<std::array<double, 3>> & /*one_over_band_width*/,
                           const bool /*domain_decomposed*/,
                           const double /*discontinuity_cutoff*/) const {

  Insist(false, "kde::reconstruction has not been implemented for this coordinate system and or "
                "dimension combination");
  return std::vector<double>(1, 0.0);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief
 * epan_kernel basis function used during reconstruction
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

//------------------------------------------------------------------------------------------------//
/*!
 * \brief
 * log_transform pushes the distribution to log space which is less sensitive to sharp derivatives.
 *
 *
 * \param[in] value of original distribution
 * \param[in] bias used to ensure positivity
 *
 * \return the logarithmic transform of the original value
 *
 * Test of kde.
 */
template <int coord> double kde<coord>::log_transform(const double value, const double bias) const {
  Require(value + bias > 0.0);
  return log(value + bias);
}

//------------------------------------------------------------------------------------------------//
/*!
 * \brief
 * log_inv_transform bring the distribution back to the original space.
 *
 *
 * \param[in] log_value of original distribution
 * \param[in] bias used to ensure positivity
 *
 * \return the logarithmic transform of the original value
 *
 * Test of kde.
 */

template <int coord>
double kde<coord>::log_inv_transform(const double log_value, const double bias) const {
  Require(log_value <= 0.0);
  return exp(log_value) - bias;
}

} // end namespace  rtt_kde

#endif // kde_kde_i_hh

//------------------------------------------------------------------------------------------------//
// end of <pkg>/kde.i.hh
//------------------------------------------------------------------------------------------------//
