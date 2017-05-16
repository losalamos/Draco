//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/Compton.cc
 * \author Kendra Keady
 * \date   Tues Feb 21 2017
 * \brief  Implementation file for compton CSK_generator interface
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

// headers provided in draco:
#include "Compton.hh"
#include "ds++/Assert.hh"
// headers provided in Compton include directory:
#include "compton_file.hh"
#include "multigroup_data_types.hh"
#include "multigroup_lib_builder.hh"

namespace rtt_compton {

// ------------ //
// Constructors //
// ------------ //

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for an existing multigroup libfile.
 *
 * This calls CSK_generator methods to read the data file and store everything in a
 * Compton data object, a smart pointer to which is then passed to (and held by)
 * the CSK_generator etemp_interp class.
 *
 * \param filehandle The name of the multigroup file to use for Compton scatters
 */
Compton::Compton(const std::string &filehandle) {

  // Check input validity
  Require(std::ifstream(filehandle).good());

  // Make a compton file object to read the multigroup data
  compton_file Cfile(false);

  // initialize the electron temperature interpolator with the mg compton data
  ei.reset(new etemp_interp(Cfile.read_mg_csk_data(filehandle)));

  // Make sure the SP exists...
  Ensure(ei);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for an existing pointwise file and a multigroup structure.
 *
 * This calls CSK_generator methods to read the pointwise library and construct a
 * multigroup Compton data object, a smart pointer to which is then passed
 * to (and held by) the CSK_generator etemp_interp class.
 *
 * \param filehandle The name of the pointwise lib to build MG data from
 * \param grp_bds    A vector containing the multigroup bounds (in keV)
 * \param n_xi       The number of angular points/Legendre moments desired
 */
Compton::Compton(const std::string &filehandle,
                 const std::vector<double> &grp_bds, const size_t nxi) {

  // Check input validity
  Require(std::ifstream(filehandle).good());
  Require(grp_bds.size() > 0);

  std::cout << "*********************************************************\n"
            << "WARNING! Building a multigroup library from scratch might\n"
            << " take a LOOOOOOONG time! (Don't say I didn't warn you.)  \n"
            << "*********************************************************\n"
            << std::endl;

  // make a group_data struct to pass to the lib builder:
  // TODO: How do we actually want to handle the weighting function?
  // Allow the user to pass it in? Make some intelligent decision at runtime?
  multigroup::Group_data grp_data = {multigroup::Library_type::EXISTING,
                                     multigroup::Weighting_function::PLANCK,
                                     filehandle, nxi, grp_bds};

  // Construct a multigroup library builder:
  multigroup_lib_builder MG_builder(grp_data);

  // build the library:
  MG_builder.build_library();

  // initialize the electron temperature interpolator with the mg compton data
  ei.reset(new etemp_interp(MG_builder.package_data()));

  // Make sure the SP exists...
  Ensure(ei);
}

// ------------ //
//  Interfaces  //
// ------------ //

//---------------------------------------------------------------------------//
/*!
 * \brief Interpolated data to a given SCALED electron temperature (T / m_e)
 *
 * This method interpolates MG Compton lib data to a given electron temperature.
 * It returns the interpolated values for ALL g, g', and angular points
 * in the specified multigroup structure.
 *
 * \param etemp The SCALED electron temperature ( temp / electron rest-mass )
 * \return      n_grp x n_grp x n_xi interpolated scattering kernel values
 */
std::vector<std::vector<std::vector<double>>>
Compton::interpolate(const double etemp) {

  // Be sure the passed electron temperature is within the bounds of the lib!
  Require(etemp >= ei->get_min_etemp());
  Require(etemp <= ei->get_max_etemp());

  // call the appropriate routine in the electron interp object
  return ei->interpolate_etemp(etemp);
}

std::vector<std::vector<std::vector<double>>>
Compton::interpolate(const double etemp) const {

  // Be sure the passed electron temperature is within the bounds of the lib!
  Require(etemp >= ei->get_min_etemp());
  Require(etemp <= ei->get_max_etemp());

  // call the appropriate routine in the electron interp object
  return ei->interpolate_etemp(etemp);
}
}
