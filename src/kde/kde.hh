//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.hh
 * \author Mathew Cleveland
 * \brief  Define class kernel density estimator class
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

// clang-format off

#ifndef kde_kde_hh
#define kde_kde_hh

#include <vector>
#include <array>
#include "c4/global.hh"

namespace kde {

//================================================================================================//
/*!
 * \class kde - kernel density estimator class for generated smoothed
 * reconstructions of point wise PDF data
 * \brief
 *
 * Returns a KDE reconstruction of a multidimensional distribution
 * 
 */
//================================================================================================//

  enum kde_coordinates{CART, CYL, SPH};

  template <int coord> class kde {
  public:
    // NESTED CLASSES AND TYPEDEFS

    // CREATORS

    // ACCESSORS

    // SERVICES
    
    //! Reconstruct distribution
    template <int dim=1>
    std::vector<double> reconstruction(const std::vector<double> &distribution, const std::vector<std::array<double,3>> &position, const std::vector<std::array<double,3>> &band_width, const bool domain_decomposed) const;

    //! Mean reconstruction 
    template <int dim=1>
    std::vector<double> mean_reconstruction(const std::vector<double> &distribution, const std::vector<std::array<double,3>> &position, const std::vector<std::array<double,3>> &band_width, const bool domain_decomposed) const;

    // STATICS

    //! Epanechikov Kernel
    double epan_kernel(const double x) const;

  protected:
    // IMPLEMENTATION

  private:
    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA
  };

template<>
template<>
std::vector<double> kde<kde_coordinates::CART>::reconstruction<1>(const std::vector<double> &distribution, const std::vector<std::array<double,3>> &position, const std::vector<std::array<double,3>> &band_width, const bool dd) const;

template<>
template<>
std::vector<double> kde<kde_coordinates::CART>::mean_reconstruction<1>(const std::vector<double> &distribution, const std::vector<std::array<double,3>> &position, const std::vector<std::array<double,3>> &band_width, const bool dd) const;



} // end namespace  kde

#include "kde.i.hh"

#endif // kde_kde_hh

//------------------------------------------------------------------------------------------------//
// end of kde/kde.hh
//------------------------------------------------------------------------------------------------//

