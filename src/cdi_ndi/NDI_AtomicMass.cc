//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/NDI_AtomicMass.cc
 * \author Ben R. Ryan
 * \date   2020 Mar 6
 * \brief  NDI_AtomicMass class declaration.
 * \note   Copyright (C) 2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "NDI_AtomicMass.hh"
#include "ndi.h"
#include "cdi_ndi/config.h"
#include "ds++/Assert.hh"
#include "ds++/path.hh"
#include <array>

namespace rtt_cdi_ndi {

//----------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------//
/*!
 * \brief Constructor for NDI atomic mass weight reader, using default path to
 *        NDI gendir file.
 */
NDI_AtomicMass::NDI_AtomicMass()
    : gendir_path(rtt_dsxx::getFilenameComponent(
          std::string(NDI_ROOT_DIR) + "/share/gendir.all",
          rtt_dsxx::FilenameComponent::FC_NATIVE)) {
  Require(rtt_dsxx::fileExists(gendir_path));
}

/*!
 * \brief Constructor for NDI atomic mass weight reader, using custom path to
 *        NDI gendir file.
 * \param[in] gendir_path_in path to gendir file
 */
NDI_AtomicMass::NDI_AtomicMass(const std::string &gendir_path_in)
    : gendir_path(gendir_path_in) {
  Require(rtt_dsxx::fileExists(gendir_path));
}

//----------------------------------------------------------------------------//
/*!
 * \brief Get atomic mass weight of an isotope with given ZAID. Use method due
 *        to T. Saller that invokes multigroup_neutron dataset which includes
 *        atomic weights.
 * \param[in] zaid ZAID of isotope for which to return the atomic mass.
 * \return mass of isotope in grams
 */
double NDI_AtomicMass::get_amw(const int zaid) const {
  // The below NDI method doesn't work for electrons
  if (zaid == -1) {
    return pc.electronMass();
  }

  int gendir_handle = -1;
  int ndi_error = -9999;
  ndi_error = NDI2_open_gendir(&gendir_handle, gendir_path.c_str());
  Require(ndi_error == 0);
  Insist(gendir_handle != -1, "gendir_handle still has default value!");

  ndi_error = NDI2_set_option_gendir(gendir_handle, NDI_LIB_TYPE_DEFAULT,
                                     "multigroup_neutron");
  Require(ndi_error == 0);

  ndi_error =
      NDI2_set_option_gendir(gendir_handle, NDI_LIBRARY_DEFAULT, "mendf71x");
  Require(ndi_error == 0);

  std::string zaid_formatted = std::to_string(zaid) + ".";

  int size = NDI2_get_size_x(gendir_handle, NDI_AT_WGT, zaid_formatted.c_str(),
                             &ndi_error);
  Require(ndi_error == 0);
  Insist(size == 1, "NDI returned more or fewer than one atomic weight?");

  std::array<double, 1> arr;
  ndi_error = NDI2_get_float64_vec_x(gendir_handle, NDI_AT_WGT,
                                     zaid_formatted.c_str(), arr.data(), size);
  Require(ndi_error == 0);

  ndi_error = NDI2_close_gendir(gendir_handle);
  Require(ndi_error == 0);

  return arr[0] * pc.amu();
}

} // namespace rtt_cdi_ndi

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_AtomicMass.cc
//----------------------------------------------------------------------------//
