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
#include "ds++/DracoStrings.hh"
#include "ds++/Query_Env.hh"
#include "ds++/SystemCall.hh"
#include "ds++/dbc.hh"

namespace rtt_cdi_ndi {

//----------------------------------------------------------------------------//
// Auxillary functions
//----------------------------------------------------------------------------//

#ifdef NDI_FOUND
/*!
 * \brief Warn if NDI library version doesn't match GENDIR version to 2 digits.
 *
 * \param[in] gendir path to the gendir file
 *
 * This function is a no-op if \c NDI_FOUND is false.
 */
void NDI_Base::warn_ndi_version_mismatch(std::string const &gendir) {

  std::string gendir_ver = rtt_dsxx::extract_version(gendir, 2);
  std::string ndi_ver = rtt_dsxx::extract_version(NDI_ROOT_DIR, 2);
  if (gendir_ver.size() == 0)
    gendir_ver = "version unknown";
  if (ndi_ver.size() == 0)
    ndi_ver = "version unknown";
  if (gendir_ver != ndi_ver) {
    using DT = Term::DracoTerminal;
    std::cout << "\n"
              << Term::ccolor(DT::error) << "WARNING: In the cdi_ndi/NDI_Base "
              << "constructor, the NDI library version (" << ndi_ver << ") is "
              << "different than the NDI GENDIR version (" << gendir_ver
              << "). \n"
              << Term::ccolor(DT::reset) << std::endl;
  }
}
#else
/*!
 * \brief Warn if NDI library version doesn't match GENDIR version to 2 digits.
 *        No-op when NDI_FOUND is false. */
void NDI_Base::warn_ndi_version_mismatch(std::string const & /*gendir*/) {}
#endif

//----------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------//

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
  warn_ndi_version_mismatch(gendir);

  for (size_t i = 0; i < mg_e_bounds.size(); i++) {
    mg_e_bounds[i] /= 1000.; // keV -> MeV
  }

  // Check that mg_e_bounds is monotonically decreasing (NDI requirement)
  for (size_t i = 1; i < mg_e_bounds.size(); i++) {
    Insist(mg_e_bounds[i] < mg_e_bounds[i - 1], "Non-monotonic mg bounds!");
  }
  Insist(mg_e_bounds[mg_e_bounds.size() - 1] > 0, "Negative mg bounds!");
}

//============================================================================//
// Stubbed implementation when NDI is unavailable
//============================================================================//

#ifndef NDI_FOUND

//! Constructor for generic NDI reader- throws when NDI not available
NDI_Base::NDI_Base(const std::string & /*dataset_in*/,
                   const std::string & /*library_in*/,
                   const std::string & /*reaction_in*/,
                   const std::vector<double> /*mg_e_bounds_in*/) {
  Insist(0, "NDI default gendir path only available when NDI is found.");
}

#else

//============================================================================//
// Normal implementation
//============================================================================//

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
    : gendir(rtt_dsxx::get_env_val<std::string>("NDI_GENDIR_PATH").second),
      dataset(dataset_in), library(library_in), reaction(reaction_in),
      mg_e_bounds(mg_e_bounds_in) {

  Require(rtt_dsxx::fileExists(gendir));

  Require(gendir.length() > 0);
  Require(dataset.length() > 0);
  Require(library.length() > 0);
  Require(reaction.length() > 0);
  Require(mg_e_bounds.size() > 0);
  warn_ndi_version_mismatch(gendir);

  for (size_t i = 0; i < mg_e_bounds.size(); i++) {
    mg_e_bounds[i] /= 1000.; // keV -> MeV
  }

  // Check that mg_e_bounds is monotonically decreasing (NDI requirement)
  Require(rtt_dsxx::is_strict_monotonic_decreasing(mg_e_bounds.begin(),
                                                   mg_e_bounds.end()));
  Require(mg_e_bounds.back() > 0);
}

#endif

} // namespace rtt_cdi_ndi

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_Base.cc
//----------------------------------------------------------------------------//
