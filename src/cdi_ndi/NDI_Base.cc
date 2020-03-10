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
#include "ds++/dbc.hh"

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
    : gendir(rtt_dsxx::getFilenameComponent(
          std::string(NDI_ROOT_DIR) + "share/gendir.all",
          rtt_dsxx::FilenameComponent::FC_NATIVE)),
      dataset(dataset_in), library(library_in), reaction(reaction_in),
      mg_e_bounds(mg_e_bounds_in) {

  Require(rtt_dsxx::fileExists(gendir));

  Require(gendir.length() > 0);
  Require(dataset.length() > 0);
  Require(library.length() > 0);
  Require(reaction.length() > 0);
  Require(mg_e_bounds.size() > 0);

  for (size_t i = 0; i < mg_e_bounds.size(); i++) {
    mg_e_bounds[i] /= 1000.; // keV -> MeV
  }

  // Check that mg_e_bounds is monotonically decreasing (NDI requirement)
  Require(rtt_dsxx::is_strict_monotonic_decreasing(mg_e_bounds.begin(),
                                                   mg_e_bounds.end()));
  Require(mg_e_bounds.back() > 0);
}

/*!
 * \brief Constructor for generic NDI reader, to be inherited by readers for
 *        specific gendir file path and dataset.
 *
 * This base constructor only sets some data members based on constructor input.
 * For more details on NDI, see
 * https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html
 *
 * \param[in] gendir_in path to non-standard NDI gendir file
 * \param[in] dataset_in name of requested dataset (provided by inherited class)
 * \param[in] library_in name of requested NDI data library
 * \param[in] reaction_in name of requested reaction
 * \param[in] mg_e_bounds_in multigroup energy bin boundaries (keV)
 */
NDI_Base::NDI_Base(const std::string &gendir_in, const std::string &dataset_in,
                   const std::string &library_in,
                   const std::string &reaction_in,
                   const std::vector<double> mg_e_bounds_in)
    : gendir(gendir_in), dataset(dataset_in), library(library_in),
      reaction(reaction_in), mg_e_bounds(mg_e_bounds_in) {

  Require(rtt_dsxx::fileExists(gendir));

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
    Insist(mg_e_bounds[i] < mg_e_bounds[i - 1], "Non-monotonic mg bounds!");
  }
  Insist(mg_e_bounds[mg_e_bounds.size() - 1] > 0, "Negative mg bounds!");
}

} // namespace rtt_cdi_ndi

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_Base.cc
//----------------------------------------------------------------------------//
