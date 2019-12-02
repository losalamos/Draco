//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyOdfmgOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:17:16 2001
 * \brief  DummyOdfmgOpacity templated class implementation file.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "DummyOdfmgOpacity.hh"
#include "ds++/Assert.hh"
#include <cmath> // pow(x,n)

namespace rtt_cdi_test {

// ------------ //
// Constructors //
// ------------ //

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for DummyOdfmgOpacity object.
 *
 * \sa The constructor assigns fixed values for all of the member
 *     data.  Every instance of this object has the same member
 *     data.
 * \code
 *     Temperatures     = { 1.0, 2.0, 3.0 }
 *     Densities        = { 0.1, 0.2 }
 *     EnergyBoundaries = { 0.05, 0.5, 5.0, 50.0 }
 *     BandBoundaries   = { 0.00, 0.125, 0.25, 0.50, 1.00 }
 * \endcode
 */
DummyOdfmgOpacity::DummyOdfmgOpacity(rtt_cdi::Reaction reaction,
                                     rtt_cdi::Model model)
    : dataFilename("none"), dataDescriptor("DummyOdfmgOpacity"),
      energyPolicyDescriptor("Odfmg"), numTemperatures(3), numDensities(2),
      numGroupBoundaries(4), numBandBoundaries(5), groupBoundaries(),
      bandBoundaries(), temperatureGrid(), densityGrid(),
      reaction_type(reaction), model_type(model) {
  temperatureGrid.resize(numTemperatures);
  densityGrid.resize(numDensities);
  groupBoundaries.resize(numGroupBoundaries);
  bandBoundaries.resize(numBandBoundaries);

  for (size_t i = 0; i < numTemperatures; ++i)
    temperatureGrid[i] = static_cast<double>(i + 1);
  for (size_t i = 0; i < numDensities; ++i)
    densityGrid[i] = static_cast<double>(i + 1) * 0.1;

  // due to the way the CDI test is implemented, these group boundaries MUST
  // match the multigroup frequency boundaries
  for (size_t i = 0; i < numGroupBoundaries; ++i)
    groupBoundaries[i] = 5.0 * std::pow(10.0, static_cast<double>(i) - 2.0);

  bandBoundaries[0] = 0.0;
  for (size_t i = 1; i < numBandBoundaries; ++i)
    bandBoundaries[i] = std::pow(2.0, static_cast<double>(i) - 4.0);
}

//---------------------------------------------------------------------------//
// Constructor for entering a different group boundary structure than the
// default
DummyOdfmgOpacity::DummyOdfmgOpacity(rtt_cdi::Reaction reaction,
                                     rtt_cdi::Model model,
                                     size_t num_groupboundaries,
                                     size_t num_bandboundaries)
    : dataFilename("none"), dataDescriptor("DummyOdfmgOpacity"),
      energyPolicyDescriptor("Odfmg"), numTemperatures(3), numDensities(2),
      numGroupBoundaries(num_groupboundaries),
      numBandBoundaries(num_bandboundaries), groupBoundaries(),
      bandBoundaries(), temperatureGrid(), densityGrid(),
      reaction_type(reaction), model_type(model) {
  temperatureGrid.resize(numTemperatures);
  densityGrid.resize(numDensities);
  groupBoundaries.resize(numGroupBoundaries);
  bandBoundaries.resize(numBandBoundaries);

  for (size_t i = 0; i < numTemperatures; ++i)
    temperatureGrid[i] = static_cast<double>(i + 1);
  for (size_t i = 0; i < numDensities; ++i)
    densityGrid[i] = static_cast<double>(i + 1) * 0.1;
  for (size_t i = 0; i < numGroupBoundaries; ++i)
    groupBoundaries[i] = 5.0 * std::pow(10.0, static_cast<double>(i) - 2.0);

  bandBoundaries[0] = 0.0;
  for (size_t i = 1; i < numBandBoundaries; ++i)
    bandBoundaries[i] =
        std::pow(2.0, static_cast<double>(i) -
                          static_cast<double>(numBandBoundaries) - 1.0);
  return;
}

// --------- //
// Accessors //
// --------- //

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of opacities (one
 *     for each group) that corresponds to the provided
 *     temperature and density.
 *
 * \code
 *     Opacity = 2 * ( temperature + density/1000 )
 *                 / ( E_high + E_low ) * 10^(band - 2)
 * \endcode
 */
std::vector<std::vector<double>>
DummyOdfmgOpacity::getOpacity(double targetTemperature,
                              double targetDensity) const {
  const size_t numGroups = getNumGroups();
  const size_t numBands = getNumBands();

  std::vector<std::vector<double>> opacity(numGroups);
  for (size_t group = 0; group < numGroups; group++) {
    opacity[group].resize(numBands);

    for (size_t band = 0; band < numBands; band++) {
      opacity[group][band] =
          2.0 * (targetTemperature + targetDensity / 1000.0) /
          (groupBoundaries[group] + groupBoundaries[group + 1]) *
          pow(10.0, static_cast<int>(band) - 2);
    }
  }
  return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of multigroup
 *     opacities corresponding to the provided vector of
 *     temperatures and a single density.  Each multigroup opacity
 *     is in itself a vector of numGroups opacities.
 *
 * \code
 *     Opacity = 2 * ( temperature + density/1000 )
 *                 / ( E_high + E_low )
 * \endcode
 */
std::vector<std::vector<std::vector<double>>>
DummyOdfmgOpacity::getOpacity(const std::vector<double> &targetTemperature,
                              double targetDensity) const {
  std::vector<std::vector<std::vector<double>>> opacity(
      targetTemperature.size());

  for (size_t i = 0; i < targetTemperature.size(); ++i) {
    opacity[i] = getOpacity(targetTemperature[i], targetDensity);
  }
  return opacity;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Opacity accessor that returns a vector of multigroup
 *     opacities corresponding to the provided vector of
 *     densities and a single temperature.  Each multigroup opacity
 *     is in itself a vector of numGroups opacities.
 *
 * \code
 *     Opacity = 2 * ( temperature + density/1000 )
 *                 / ( E_high + E_low )
 * \endcode
 */
std::vector<std::vector<std::vector<double>>>
DummyOdfmgOpacity::getOpacity(double targetTemperature,
                              const std::vector<double> &targetDensity) const {
  std::vector<std::vector<std::vector<double>>> opacity(targetDensity.size());

  //call our regular getOpacity function for every target density
  for (size_t i = 0; i < targetDensity.size(); ++i) {
    opacity[i] = getOpacity(targetTemperature, targetDensity[i]);
  }
  return opacity;
}

} // end namespace rtt_cdi_test

//---------------------------------------------------------------------------//
// end of DummyOdfmgOpacity.cc
//---------------------------------------------------------------------------//
