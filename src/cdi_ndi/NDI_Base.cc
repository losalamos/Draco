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
#include "ds++/SystemCall.hh"
#include "ds++/dbc.hh"

namespace rtt_cdi_ndi {

//----------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------//

//============================================================================//
// Stubbed implementation when NDI is unavailable
//============================================================================//

#ifndef NDI_FOUND

//! Constructor for generic NDI reader- throws when NDI not available
NDI_Base::NDI_Base(const std::string & /*dataset_in*/,
                   const std::string & /*library_in*/) {
  Insist(0, "NDI default gendir path only available when NDI is found.");
}

#else

//============================================================================//
// Normal implementation
//============================================================================//

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
 */
NDI_Base::NDI_Base(const std::string &gendir_in, const std::string &dataset_in,
                   const std::string &library_in)
    : gendir(gendir_in), dataset(dataset_in), library(library_in) {

  Require(rtt_dsxx::fileExists(gendir));

  Require(gendir.length() > 0);
  Require(dataset.length() > 0);
  Require(library.length() > 0);
}

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
 */
NDI_Base::NDI_Base(const std::string &dataset_in, const std::string &library_in)
    : gendir(rtt_dsxx::getFilenameComponent(
          std::string(NDI_DATA_DIR) + rtt_dsxx::dirSep + "gendir",
          rtt_dsxx::FilenameComponent::FC_NATIVE)),
      dataset(dataset_in), library(library_in) {

  Require(rtt_dsxx::fileExists(gendir));

  Require(gendir.length() > 0);
  Require(dataset.length() > 0);
  Require(library.length() > 0);
}

#endif

} // namespace rtt_cdi_ndi

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_Base.cc
//----------------------------------------------------------------------------//
