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

  //! one-line description
  template <int coord>
  template <int dim>
  void kde<coord>::reconstruction(const std::vector<double> & /*distribution*/,
                                     const std::vector<double> /*position*/,
                                     const std::vector<double> /*band_width*/,
                                     std::vector<double> & /*final_distribution*/) const {

  Insist(false, "kde::reconstruction has not been implemented for this coordinate system and or "
                "dimension combination");
  }


} // end namespace  kde

#endif // kde_kde_i_hh

//------------------------------------------------------------------------------------------------//
// end of <pkg>/kde.i.hh
//------------------------------------------------------------------------------------------------//

