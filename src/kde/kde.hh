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
    void reconstruction(const std::vector<double> &distribution, const std::vector<double> position, const std::vector<double> band_width, std::vector<double> &final_distribution) const;

    // STATICS

  protected:
    // IMPLEMENTATION

  private:
    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    // DATA
  };

template<>
template<>
void kde<kde_coordinates::CART>::reconstruction<1>(const std::vector<double> &distribution, const std::vector<double> position, const std::vector<double> band_width, std::vector<double> &final_distribution) const;
template<>
template<>
void kde<kde_coordinates::CART>::reconstruction<2>(const std::vector<double> &distribution, const std::vector<double> position, const std::vector<double> band_width, std::vector<double> &final_distribution) const;


} // end namespace  kde

#include "kde.i.hh"

#endif // kde_kde_hh

//------------------------------------------------------------------------------------------------//
// end of kde/kde.hh
//------------------------------------------------------------------------------------------------//

