//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   kde/kde.i.hh
 * \author Mathew Cleveland
 * \date   Nov. 10th 2020
 * \brief  Member definitions of class kde
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

// clang-format off

#ifndef kde_kde_i_hh
#define kde_kde_i_hh

#include "ds++/dbc.hh"

namespace kde {

  //! DEFAULT constructor to return error if instantiation is not found
  template <int coord>
  template <int dim>
  std::vector<double> kde<coord>::reconstruction(const std::vector<double> & /*distribution*/,
                                     const std::vector<std::array<double,3>> &/*position*/,
                                     const std::vector<std::array<double,3>> &/*band_width*/,
                                     const bool /*domain_decomposed*/) const {

  Insist(false, "kde::reconstruction has not been implemented for this coordinate system and or "
                "dimension combination");
  return std::vector<double>(1,0.0);
  }

  //! DEFAULT constructor to return error if instantiation is not found
  template <int coord>
  template <int dim>
  std::vector<double> kde<coord>::mean_reconstruction(const std::vector<double> & /*distribution*/,
                                     const std::vector<std::array<double,3>> &/*position*/,
                                     const std::vector<std::array<double,3>> &/*band_width*/,
                                     const bool /*domain_decomposed*/) const {

  Insist(false, "kde::reconstruction has not been implemented for this coordinate system and or "
                "dimension combination");
  return std::vector<double>(1,0.0);
  }



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
  template<int coord>
  double kde<coord>::epan_kernel(const double x) const{
      const double x2 = x*x;
      const double result = 0.75*(1.0-x2);
      if (x2>1.0)
          return 0.0;
      else
        return result;
  }

} // end namespace  kde

#endif // kde_kde_i_hh

//------------------------------------------------------------------------------------------------//
// end of <pkg>/kde.i.hh
//------------------------------------------------------------------------------------------------//

