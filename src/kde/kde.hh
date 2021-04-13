//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.hh
 * \author Mathew Cleveland
 * \brief  Define class kernel density estimator class
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef kde_kde_hh
#define kde_kde_hh

#include "quick_index.hh"
#include "c4/global.hh"
#include <array>
#include <vector>

namespace rtt_kde {

enum kde_coordinates { CART, CYL, SPH };

//================================================================================================//
/*!
 * \class kde
 * \brief kernel density estimator class for generated smoothed reconstructions of point wise PDF
 *        data
 * 
 * \tparam coord enumeration specifying the KDE coordinate system to use.
 *
 * Returns a KDE reconstruction of a multidimensional distribution
 */
//================================================================================================//
template <int coord> class kde {
public:
  // NESTED CLASSES AND TYPEDEFS

  // CREATORS

  // ACCESSORS

  // SERVICES

  //! Reconstruct distribution
  template <int dim = 1>
  std::vector<double> reconstruction(const std::vector<double> &distribution,
                                     const std::vector<std::array<double, 3>> &one_over_band_width,
                                     const quick_index<dim> &qindex,
                                     const double discontinuity_cutoff = 0.0) const;

  //! Reconstruct distribution in logarithmic space
  template <int dim = 1>
  std::vector<double>
  log_reconstruction(const std::vector<double> &distribution,
                     const std::vector<std::array<double, 3>> &one_over_band_width,
                     const quick_index<dim> &qindex, const double discontinuity_cutoff = 0.0) const;

  // STATICS

  //! Epanechikov Kernel
  double epan_kernel(const double x) const;

  //! Transform the solution into log space
  double log_transform(const double value, const double bias) const;

  //! Move the solution back from log space
  double log_inv_transform(const double log_value, const double bias) const;

protected:
  // IMPLEMENTATION

private:
  // NESTED CLASSES AND TYPEDEFS

  // IMPLEMENTATION

  // DATA
};

//! Forward declaration of the reconstruction 1D Cartesian reconstruction.
template <>
template <>
std::vector<double> kde<kde_coordinates::CART>::reconstruction<1>(
    const std::vector<double> &distribution,
    const std::vector<std::array<double, 3>> &one_over_band_width, const quick_index<1> &qindex,
    const double discontinuity_cutoff) const;

template <>
template <>
std::vector<double> kde<kde_coordinates::CART>::log_reconstruction<1>(
    const std::vector<double> &distribution,
    const std::vector<std::array<double, 3>> &one_over_band_width, const quick_index<1> &qindex,
    const double discontinuity_cutoff) const;

} // end namespace rtt_kde

#include "kde.i.hh"

#endif // kde_kde_hh

//------------------------------------------------------------------------------------------------//
// end of kde/kde.hh
//------------------------------------------------------------------------------------------------//
