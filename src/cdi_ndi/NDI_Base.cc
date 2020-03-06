//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/NDI_Base.cc
 * \author Ben R. Ryan
 * \date   2020 Feb 4
 * \brief  NDI_Base member definitions.
 * \note   Copyright (C) 2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "NDI_Base.hh"

namespace rtt_cdi_ndi {

//----------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------//
/*!
 * \brief Constructor for generic NDI reader, to be inherited by readers for
 *        specific dataset.
 *
 * This base constructor only sets some data members based on constructor input.
 * For more details on NDI, see
 * https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html
 *
 * \param[in] dataset_in name of requested dataset (provided by inherited class)
 * \param[in] library_in name of requested NDI data library
 * \param[in] reaction_in name of requested reaction
 * \param[in] mg_e_bounds_in multigroup energy bin boundaries (keV)
 */
NDI_Base::NDI_Base(const std::string &dataset_in, const std::string &library_in,
                   const std::string &reaction_in,
                   const std::vector<double> mg_e_bounds_in)
    : dataset(dataset_in), library(library_in), reaction(reaction_in),
      mg_e_bounds(mg_e_bounds_in) {
  gendir = std::string(NDI_ROOT_DIR) + "share/gendir.all";

  Require(gendir.length() > 0);
  Require(dataset.length() > 0);
  Require(library.length() > 0);
  Require(reaction.length() > 0);
  Require(mg_e_bounds.size() > 0);

  for (size_t i = 0; i < mg_e_bounds.size(); i++) {
    mg_e_bounds[i] /= 1000.; // keV -> MeV
  }

  // Check that mg_e_bounds is monotonically decreasing (NDI requirement)
  for (size_t i = 1; i < mg_e_bounds.size(); i++) {
    Require(mg_e_bounds[i] < mg_e_bounds[i - 1]);
  }
  Require(mg_e_bounds[mg_e_bounds.size() - 1] > 0);
}

} // namespace rtt_cdi_ndi

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_Base.cc
//----------------------------------------------------------------------------//
