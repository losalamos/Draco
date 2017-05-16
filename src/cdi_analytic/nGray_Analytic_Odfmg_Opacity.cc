//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/nGray_Analytic_Odfmg_Opacity.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 13 11:19:59 2001
 * \brief  nGray_Analytic_Odfmg_Opacity class member definitions.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "nGray_Analytic_Odfmg_Opacity.hh"
#include "nGray_Analytic_MultigroupOpacity.hh"
#include "ds++/Packing_Utils.hh"
#include "ds++/dbc.hh"

namespace rtt_cdi_analytic {

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic multigroup opacity model.
 *
 * This constructor builds an opacity model defined by the
 * rtt_cdi_analytic::nGray_Analytic_Opacity_Model derived class argument.
 *
 * The reaction type for this instance of the class is determined by the
 * rtt_cdi::Reaction argument.
 *
 * The group structure (in keV) must be provided by the groups argument.  The
 * number of nGray_Analytic_Opacity_Model objects given in the models argument
 * must be equal to the number of groups.
 *
 * \param groups vector containing the group boundaries in keV from lowest to
 * highest
 *
 * \param models vector containing shared_ptrs to Analytic_Model derived types
 * for each group, the size should be groups.size() - 1
 *
 * \param reaction_in rtt_cdi::Reaction type (enumeration)
 *
 */
nGray_Analytic_Odfmg_Opacity::nGray_Analytic_Odfmg_Opacity(
    const sf_double &groups, const sf_double &bands,
    const sf_Analytic_Model &models, rtt_cdi::Reaction reaction_in,
    rtt_cdi::Model model_in)
    : Analytic_Odfmg_Opacity(groups, bands, reaction_in, model_in),
      group_models(models) {
  Require(models.size() == groups.size() - 1);
  Require(
      rtt_dsxx::is_strict_monotonic_increasing(groups.begin(), groups.end()));
  Require(rtt_dsxx::is_strict_monotonic_increasing(bands.begin(), bands.end()));
}

//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 *
 * This constructor rebuilds and nGray_Analytic_Odfmg_Opacity from a
 * vector<char> that was created by a call to pack().  It can only rebuild
 * Analytic_Model types that have been registered in the
 * rtt_cdi_analytic::Opacity_Models enumeration.
 */
nGray_Analytic_Odfmg_Opacity::nGray_Analytic_Odfmg_Opacity(
    const sf_char &packed)
    : Analytic_Odfmg_Opacity(packed), group_models() {
  // the packed size must be at least 5 integers (number of groups, number of
  // bands, reaction type, model type, analytic model indicator)
  Require(packed.size() >= 5 * sizeof(int));

  unsigned const base_size = Analytic_Odfmg_Opacity::packed_size();

  // make an unpacker
  rtt_dsxx::Unpacker unpacker;

  // register the unpacker
  unpacker.set_buffer(packed.size() - base_size, &packed[base_size]);

  // unpack the number of group boundaries
  sf_double const &group_boundaries = getGroupBoundaries();
  int ngrp_bounds = group_boundaries.size();
  int num_groups = ngrp_bounds - 1;

  // make the group boundaries and model vectors
  group_models.resize(num_groups);

  // now unpack the models
  std::vector<sf_char> models(num_groups);
  int model_size = 0;
  for (size_t i = 0; i < models.size(); i++) {
    // unpack the size of the analytic model
    unpacker >> model_size;
    Check(static_cast<size_t>(model_size) >= sizeof(int));

    models[i].resize(model_size);

    // unpack the model
    for (size_t j = 0; j < models[i].size(); j++)
      unpacker >> models[i][j];
  }

  // now rebuild the analytic models
  int indicator = 0;
  for (size_t i = 0; i < models.size(); i++) {
    // reset the buffer
    unpacker.set_buffer(models[i].size(), &models[i][0]);

    // get the indicator for this model (first packed datum)
    unpacker >> indicator;

    // now determine which analytic model we need to build
    if (indicator == rtt_cdi_analytic::CONSTANT_ANALYTIC_OPACITY_MODEL) {
      group_models[i].reset(new Constant_Analytic_Opacity_Model(models[i]));
    } else if (indicator ==
               rtt_cdi_analytic::POLYNOMIAL_ANALYTIC_OPACITY_MODEL) {
      group_models[i].reset(new Polynomial_Analytic_Opacity_Model(models[i]));
    } else if (indicator ==
               rtt_cdi_analytic::STIMULATED_EMISSION_ANALYTIC_OPACITY_MODEL) {
      group_models[i].reset(
          new Stimulated_Emission_Analytic_Opacity_Model(models[i]));
    } else {
      Insist(false, "Unregistered analytic opacity model!");
    }

    Ensure(group_models[i]);
  }
}

//---------------------------------------------------------------------------//
// OPACITY INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return the group opacities given a scalar temperature and density.
 *
 * Given a scalar temperature and density, return the group opacities
 * (vector<double>) for the reaction type specified by the constructor.  The
 * analytic opacity model is specified in the constructor
 * (nGray_Analytic_Odfmg_Opacity()).
 *
 * \param temperature material temperature in keV
 * \param density material density in g/cm^3
 * \return group opacities (coefficients) in cm^2/g
 *
 */
std::vector<std::vector<double>>
nGray_Analytic_Odfmg_Opacity::getOpacity(double targetTemperature,
                                         double targetDensity) const {
  Require(targetTemperature >= 0.0);
  Require(targetDensity >= 0.0);

  const size_t numBands = getNumBands();
  const size_t numGroups = getNumGroups();

  sf_double const &group_bounds = this->getGroupBoundaries();

  // return opacities
  std::vector<std::vector<double>> opacity(numGroups);

  // loop through groups and get opacities
  for (size_t group = 0; group < opacity.size(); group++) {
    Check(group_models[group]);

    opacity[group].resize(numBands);

    // assign the opacity based on the group model to the first band
    opacity[group][0] = group_models[group]->calculate_opacity(
        targetTemperature, targetDensity, group_bounds[group],
        group_bounds[group + 1]);

    Check(opacity[group][0] >= 0.0);

    //copy the opacity to the rest of the bands
    for (size_t band = 1; band < numBands; band++) {
      opacity[group][band] = opacity[group][0];
    }
  }

  return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of multigroupband
 *     opacity 2-D vectors that correspond to the provided vector of
 *     temperatures and a single density value.
 */
std::vector<std::vector<std::vector<double>>>
nGray_Analytic_Odfmg_Opacity::getOpacity(
    const std::vector<double> &targetTemperature, double targetDensity) const {
  std::vector<std::vector<std::vector<double>>> opacity(
      targetTemperature.size());

  for (size_t i = 0; i < targetTemperature.size(); ++i) {
    opacity[i] = getOpacity(targetTemperature[i], targetDensity);
  }
  return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of multigroupband
 *     opacity 2-D vectors that correspond to the provided
 *     temperature and a vector of density values.
 */
std::vector<std::vector<std::vector<double>>>
nGray_Analytic_Odfmg_Opacity::getOpacity(
    double targetTemperature, const std::vector<double> &targetDensity) const {
  std::vector<std::vector<std::vector<double>>> opacity(targetDensity.size());

  //call our regular getOpacity function for every target density
  for (size_t i = 0; i < targetDensity.size(); ++i) {
    opacity[i] = getOpacity(targetTemperature, targetDensity[i]);
  }
  return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack an analytic odfmg opacity.
 *
 * This function will pack up the Analytic_Mulitgroup_Opacity into a char
 * array (represented by a vector<char>).  The nGray_Analytic_Opacity_Model
 * derived class must have a pack function; this is enforced by the virtual
 * nGray_Analytic_Opacity_Model base class.
 */
nGray_Analytic_Odfmg_Opacity::sf_char
nGray_Analytic_Odfmg_Opacity::pack() const {
  // make a packer
  rtt_dsxx::Packer packer;

  // make a char array
  sf_char packed = Analytic_Odfmg_Opacity::pack();

  // first pack up models
  std::vector<sf_char> models(group_models.size());
  size_t num_bytes_models = 0;

  // loop through and pack up the models
  for (size_t i = 0; i < models.size(); i++) {
    Check(group_models[i]);

    models[i] = group_models[i]->pack();
    num_bytes_models += models[i].size();
  }

  // now add up the total size; number of groups + 1 size_t for number of
  // groups, number of bands + 1 size_t for number of
  // bands, number of models + size in each model + models, 1 size_t for
  // reaction type, 1 size_t for model type
  size_t base_size = packed.size();
  size_t size = models.size() * sizeof(int) + num_bytes_models;

  // extend the char array
  packed.resize(size + base_size);

  // set the buffer
  packer.set_buffer(size, &packed[base_size]);

  // pack each models size and data
  for (size_t i = 0; i < models.size(); i++) {
    // pack the size of this model
    packer << static_cast<int>(models[i].size());

    // now pack the model data
    for (size_t j = 0; j < models[i].size(); j++)
      packer << models[i][j];
  }

  Ensure(packer.get_ptr() == &packed[0] + size + base_size);
  return packed;
}

} // end namespace rtt_cdi_analytic

//---------------------------------------------------------------------------//
// end of nGray_Analytic_Odfmg_Opacity.cc
//---------------------------------------------------------------------------//
