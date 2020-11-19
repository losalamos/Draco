//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/NDI_CP_Eloss.cc
 * \author Ben R. Ryan
 * \date   2020 Jun 3
 * \brief  NDI_CP_Eloss member definitions.
 * \note   Copyright (C) 2019-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "NDI_CP_Eloss.hh"
#include "ds++/DracoStrings.hh"

namespace rtt_cdi_ndi {

// Protect actual NDI calls with NDI_FOUND macro:
#ifdef NDI_FOUND
//----------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------//
/*!
 * \brief Constructor for NDI reader specific to TN DEDX data with provided path
 *        to gendir file.
 *
 * \param[in] gendir_in path to gendir file
 * \param[in] library_in name of requested NDI data library
 * \param[in] reaction_in name of requested reaction
 * \param[in] mg_e_bounds_in energy boundaries of multigroup bins (keV)
 */
NDI_CP_Eloss::NDI_CP_Eloss(const std::string &gendir_in,
                               const std::string &library_in,
                               rtt_cdi::CParticle target_in,
                               rtt_cdi::CParticle projectile_in)
    : rtt_cdi::CPEloss(target_in,
      projectile_in, rtt_cdi::CPModelType::TABULAR_ETYPE,
      rtt_cdi::CPModelAngleCutoff::NONE),
      NDI_Base(gendir_in, "dedx", library_in) {
  // Set angle cutoff by parsing library_in

  load_ndi();
}
/*!
 * \brief Constructor for NDI reader specific to TN DEDX data using default
 *        gendir file.
 *
 * \param[in] library_in name of requested NDI data library
 * \param[in] reaction_in name of requested reaction
 * \param[in] mg_e_bounds_in energy boundaries of multigroup bins (keV)
 */
NDI_CP_Eloss::NDI_CP_Eloss(const std::string &library_in,
                               rtt_cdi::CParticle target_in,
                               rtt_cdi::CParticle projectile_in)
    : rtt_cdi::CPEloss(target_in,
      projectile_in, rtt_cdi::CPModelType::TABULAR_ETYPE,
      rtt_cdi::CPModelAngleCutoff::NONE),
      NDI_Base("dedx", library_in) {
  // Set angle cutoff by parsing library_in

  load_ndi();
}
//----------------------------------------------------------------------------//
/*!
 * \brief Load NDI dataset
 *
 * This function opens an NDI file, navigates to the appropriate dataset, reads
 * the data into internal buffers, and closes the file. For more details on NDI,
 * see https://xweb.lanl.gov/projects/data/nuclear/ndi/ndi.html
 */
void NDI_CP_Eloss::load_ndi() {
  int gendir_handle = -1;
  int dataset_handle = -1;
  int ndi_error = -9999;
  constexpr int c_str_len = 4096;
  char c_str_buf[c_str_len];

  // Open gendir file (index of a complete NDI dataset)
  printf("about to open %s\n", gendir.c_str());
  ndi_error = NDI2_open_gendir(&gendir_handle, gendir.c_str());
  Require(ndi_error == 0);
  Insist(gendir_handle != -1, "gendir_handle still has default value!");

  // Set dataset option by changing default value for this handle
  printf("dataset: %s\n", dataset.c_str());
  ndi_error = NDI2_set_option_gendir(gendir_handle, NDI_LIB_TYPE_DEFAULT,
                                     dataset.c_str());
  Require(ndi_error == 0);

  //! Set library option by changing default value for this handle
  printf("library: %s\n", library.c_str());
  ndi_error = NDI2_set_option_gendir(gendir_handle, NDI_LIBRARY_DEFAULT,
                                     library.c_str());
  Require(ndi_error == 0);

  //! Get dataset handle
  ndi_error = NDI2_make_handle(gendir_handle, dataset.c_str(), &dataset_handle);
  Require(ndi_error == 0);
  Insist(dataset_handle != -1, "dataset_handle still has default value!");

  //! Set projectile isotope
  ndi_error = NDI2_set_isotope(dataset_handle,
                                std::to_string(projectile.get_zaid()).c_str());
  printf("zaid: %i\n", projectile.get_zaid());

  int num_targets = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_TARGET, &num_targets);
  Require(ndi_error == 0);
  printf("num_targets: %i\n", num_targets);

  int num_grps = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_GRPS, &num_grps);
  Require(ndi_error == 0);
  printf("num_grps: %i\n", num_grps);

  printf("NDI_DEDX: %i\n", NDI_DEDX);

  int num_densities = 0;
  ndi_error = NDI2_get_int_val(dataset_handle, NDI_NUM_DENSITIES, &num_densities);
  Require(ndi_error == 0);
  printf("num_densities: %i\n", num_densities);

  printf("[%s : %i] not implemented!\n", __FILE__, __LINE__);
  exit(-1);
}
//----------------------------------------------------------------------------//
/*!
 * \brief Interpolate the tabulated stopping power for a given material and
 *        projectile state.
 * \param[in] temperature Material temperature [keV]
 * \param[in] density Material density [g cm^-3]
 * \param[in] partSpeed Particle speed [cm shk^-1]
 */
double NDI_CP_Eloss::getEloss(const double temperature,
                                  const double density,
                                  const double partSpeed) const {
                                    return 0.;
}
#endif // NDI_FOUND
} // namespace rtt_cdi_ndi

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_CP_Eloss.cc
//----------------------------------------------------------------------------//
