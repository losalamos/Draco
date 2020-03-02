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
 * \param[in] gendir_in path to gendir file
 * \param[in] dataset_in name of requested dataset (provided by inherited class)
 * \param[in] library_in name of requested NDI data library
 * \param[in] reaction_in name of requested reaction
 * \param[in] mg_form_in name of multigroup energy discretization to use
 */
NDI_Base::NDI_Base(const std::string &gendir_in, const std::string &dataset_in,
                   const std::string &library_in,
                   const std::string &reaction_in, const MG_FORM mg_form_in)
    : gendir(gendir_in), dataset(dataset_in), library(library_in),
      reaction(reaction_in), mg_form(mg_form_in) {

  Require(gendir.length() > 0);
  Require(dataset.length() > 0);
  Require(library.length() > 0);
  Require(reaction.length() > 0);

  mg_form_map.insert(mg_pair(MG_FORM::LANL4, "4_lanl"));
  mg_form_map.insert(
      mg_pair(MG_FORM::NOT_SET, "No multigroup representation set!"));
}

} // namespace rtt_cdi_ndi

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_Base.cc
//----------------------------------------------------------------------------//
