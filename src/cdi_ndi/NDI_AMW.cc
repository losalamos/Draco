//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   cdi_ndi/NDI_AMW.cc
 * \author Ben R. Ryan
 * \date   2020 Mar 6
 * \brief  NDI_AMW class declaration.
 * \note   Copyright (C) 2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "NDI_AMW.hh"
#include "ndi.h"
#include "cdi_ndi/config.h"
#include "ds++/Assert.hh"

namespace rtt_cdi_ndi {

//----------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------//
/*!
 * \brief Constructor for NDI atomic mass weight reader, using default path to
 *        NDI gendir file.
 */
NDI_AMW::NDI_AMW() {
  gendir_path = std::string(NDI_ROOT_DIR) + "/share/gendir.all";
}

/*!
 * \brief Constructor for NDI atomic mass weight reader, using custom path to
 *        NDI gendir file.
 * \param[in] gendir_path_in path to gendir file
 */
NDI_AMW::NDI_AMW(const std::string gendir_path_in) {
  gendir_path = gendir_path_in;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Get atomic mass weight of an isotope with given ZAID. Use method due
 *        to T. Saller that invokes multigroup_neutron dataset which includes
 *        atomic weights.
 * \param[in] zaid ZAID of isotope for which to return the atomic mass.
 * \return mass of isotope in grams
 */
double NDI_AMW::get_amw(const int zaid) const {
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

  int size = NDI2_get_size_x(gendir_handle, NDI_AT_WGT,
                             std::to_string(zaid).c_str(), &ndi_error);
  Require(ndi_error == 0);
  Insist(size == 1, "NDI returned more than one atomic weight?");

  double arr[1];
  ndi_error = NDI2_get_float64_vec_x(gendir_handle, NDI_AT_WGT,
                                     std::to_string(zaid).c_str(), arr, size);
  Require(ndi_error == 0);

  ndi_error = NDI2_close_gendir(gendir_handle);
  Require(ndi_error == 0);

  return arr[0] * pc.amu();
}

} // namespace rtt_cdi_ndi

//----------------------------------------------------------------------------//
// End cdi_ndi/NDI_AMW.cc
//----------------------------------------------------------------------------//
