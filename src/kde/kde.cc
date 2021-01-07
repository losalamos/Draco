//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.cc
 * \author Mathew Cleveland
 * \date   November 10th 2020
 * \brief  Explicitly defined KDE functions for various dimensions and coordinate
 *         KDE or Kernel Density Estimators are unbiased statical based
 *         reconstruction.  They can significantly increase the convergence
 *         rate of statical distributions. The KDE performs a reconstruction by
 *         evaluating a mean over some discrete kernel shape. In this DRACO
 *         implementation the mean is evaluated based on the sample locations
 *         that are bound by the kernel shape.  A renormalization is used to
 *         ensure the proper mean is returned given there is no guarantee the
 *         full kernel (which integrates exactly to 1) will be integrated fully
 *         in space. This renormalization also avoids the need for boundary
 *         fix-ups which are typically used in KDE applications to account for
 *         the kernel extending beyond the bounds of the spatial domain. Other
 *         approaches that could be considered are quadrature based approaches
 *         that fully sample the Kernel space reducing the need for the
 *         normalization.
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC., All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "kde.hh"
#include <cmath>
#include <iostream>
#include <numeric>

namespace rtt_kde {

//------------------------------------------------------------------------------------------------//
/*!
 * \brief Cartesian geometry KDE reconstruction of a 1D distribution
 * 
 * \pre The local reconstruction data is passed into this function which
 * includes the original data distribution, its spatial position, and the
 * optimal bandwidth to be used at each point.
 *
 * \param[in] distribution original data to be reconstructed
 * \param[in] position local of the original data 
 * \param[in] one_over_band_width inverse bandwidth size to be used at each data location
 * \param[in] domain_decomposed bool flag to switch between domain replicated and decomposed data
 * \return final local KDE function distribution reconstruction
 *
 * \post the local reconstruction of the original data is returned.
 */
template <>
template <>
std::vector<double> kde<kde_coordinates::CART>::reconstruction<1>(
    const std::vector<double> &distribution, const std::vector<std::array<double, 3>> &position,
    const std::vector<std::array<double, 3>> &one_over_band_width,
    const bool domain_decomposed) const {
  const int64_t local_size = distribution.size();
  Require(static_cast<int64_t>(position.size()) == local_size);
  Require(static_cast<int64_t>(one_over_band_width.size()) == local_size);

  // used for the zero accumulation conservation
  int64_t size = local_size;
  double global_conservation = std::accumulate(distribution.begin(), distribution.end(), 0.0);
  std::vector<double> result(local_size, 0.0);
  std::vector<double> normal(local_size, 0.0);
  if (domain_decomposed) {
    // minimize global values and only allocate them in DD problems
    int64_t global_lower_bound = 0;
    std::vector<double> global_distribution;
    std::vector<double> global_x_position;

    // calculate global off sets
    int n_ranks = rtt_c4::nodes();
    int64_t izero(0);
    std::vector<int64_t> rank_size(rtt_c4::nodes(), izero);
    rank_size[rtt_c4::node()] = local_size;
    rtt_c4::global_sum(rank_size.data(), n_ranks);
    size = std::accumulate(rank_size.begin(), rank_size.end(), izero);
    std::vector<int64_t> accum_rank_size(rank_size);
    std::partial_sum(rank_size.begin(), rank_size.end(), accum_rank_size.begin());

    if (rtt_c4::node() > 0) {
      global_lower_bound = accum_rank_size[rtt_c4::node() - 1];
    }

    // set up global arrays
    global_distribution.resize(size, 0.0);
    global_x_position.resize(size, 0.0);

    // build up global positions
    for (int i = 0; i < local_size; i++) {
      global_x_position[i + global_lower_bound] = position[i][0];
      global_distribution[i + global_lower_bound] = distribution[i];
    }

    rtt_c4::global_sum(global_x_position.data(), size);
    rtt_c4::global_sum(global_distribution.data(), size);
    rtt_c4::global_sum(global_conservation);

    // now apply the kernel to the local ranks
    for (int i = 0; i < local_size; i++) {
      const double x0 = position[i][0];
      const double one_over_h = one_over_band_width[i][0];
      // fetch local contribution
      for (int j = 0; j < size; j++) {
        const double x = global_x_position[j];
        const double u = (x0 - x) * one_over_h;
        const double weight = (epan_kernel(u)) * one_over_h;
        result[i] += global_distribution[j] * weight;
        normal[i] += weight;
      }
    }
  } else { // local reconstruction only
    // now apply the kernel to the local ranks
    for (int i = 0; i < local_size; i++) {
      const double x0 = position[i][0];
      const double one_over_h = one_over_band_width[i][0];
      // fetch local contribution
      for (int j = 0; j < local_size; j++) {
        const double x = position[j][0];
        const double u = (x0 - x) * one_over_h;
        const double weight = (epan_kernel(u)) * one_over_h;
        result[i] += distribution[j] * weight;
        normal[i] += weight;
      }
    }
  }

  // normalize the integrated weight contributions
  for (int i = 0; i < local_size; i++) {
    Check(normal[i] > 0.0);
    result[i] /= normal[i];
  }

  double reconstruction_conservation = std::accumulate(result.begin(), result.end(), 0.0);

  if (domain_decomposed) {
    // accumulate global contribution
    rtt_c4::global_sum(reconstruction_conservation);
  }

  if (!rtt_dsxx::soft_equiv(reconstruction_conservation, 0.0) &&
      !rtt_dsxx::soft_equiv(global_conservation, 0.0)) {
    // Totals are non-zero so scale the result for conservation
    for (int i = 0; i < local_size; i++)
      result[i] *= global_conservation / reconstruction_conservation;
  } else {
    // a zero distribution is possible. If it occurs fall back to residual conservation;
    const double res = global_conservation - reconstruction_conservation;
    for (int i = 0; i < local_size; i++)
      result[i] += res / double(size);
  }

  return result;
}

} // end namespace rtt_kde

//------------------------------------------------------------------------------------------------//
// end of kde/kde.cc
//------------------------------------------------------------------------------------------------//
