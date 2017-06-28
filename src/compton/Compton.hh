//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/Compton.hh
 * \author Kendra Keady
 * \date   Mon Feb 27 2017
 * \brief  Header file for compton CSK_generator interface
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __compton_Compton_hh__
#define __compton_Compton_hh__

// C++ standard library dependencies
#include <iostream>
#include <memory>
#include <vector>
// headers provided in compton CSK_generator include directory
#include "etemp_interp.hh"
#include "multigroup_compton_data.hh"

namespace rtt_compton {
//===========================================================================//
/*!
 * \class Compton
 *
 * \brief Provides access to relativistic Compton scattering angle and
 * multigroup frequency distributions from the CSK_generator project.
 *
 * This interface class allows the client to:
 * 1) access (interpolate) data from existing multigroup CSK_generator libraries
 * 2) build new multigroup libraries from existing CSK_generator pointwise
      libraries
 * 3) obtain auxiliary information for existing multigroup libraries
 *    (electron temperature bounds, frequency group structures, etc)
 *
 * This class is designed to be used with the CSK_generator library and headers.
 * If this is not found at the CMake configure step, the lib_compton portion of
 * draco will not be built.
 *
 * <b>User's environment</b>
 *
 * Cmake searches for the CSK_generator library/include headers during the
 * configuration step. The script that does this is located at:
 *
 * /draco/config/findCOMPTON.cmake
 *
 */

/*!
 * \example compton/test/tCompton.cc
 *
 * This unit test demonstrates the two methods for constructing a Compton
 * object, and exercises all routines for interpolation and data access.
*/

class Compton {

private:
  //! Shared pointer to an electron interpolation object:
  std::shared_ptr<etemp_interp> ei;

public:
  //! Constructor for an existing multigroup library
  explicit Compton(const std::string &filehandle);

  //! Constructor to build a multigroup library from an existing pointwise file
  Compton(const std::string &file, const std::vector<double> &group_bounds,
          const std::string &opac_type, const std::string &wt_func,
          const bool induced, const size_t n_xi = 0);

  //! Interpolation of all csk opacity data to a certain electron temperature:
  std::vector<std::vector<std::vector<std::vector<double>>>>
  interpolate_csk(const double etemp) const;

  //! Interpolation of all nu_ratio data to an electron temperature:
  std::vector<std::vector<double>>
  interpolate_nu_ratio(const double etemp) const;

  //! Retrieve group structure for the given library:
  std::vector<double> get_group_bounds() const {
    return ei->get_Cdata()->get_group_bds();
  }

  //! Retrieve min electron temperature for the given library:
  double get_min_etemp() const { return ei->get_min_etemp(); }

  //! Retrieve max electron temperature for the given library:
  double get_max_etemp() const { return ei->get_max_etemp(); }

  //! Retrieve number of groups in the given multigroup structure:
  size_t get_num_groups() const { return ei->get_Cdata()->get_n_grps(); }

  //! Retrieve number of angular moments/evaluation points in the lib data:
  size_t get_num_xi() const { return ei->get_Cdata()->get_n_xi_pts(); }

  //! Retrieve electron temperature eval points (diagnostic use)
  std::vector<double> get_etemp_pts() const { return ei->get_etemp_pts(); }
};
}

#endif
