//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_ipcress/test/cdi_ipcress_test.hh
 * \author Thomas M. Evans
 * \date   Fri Oct 12 15:36:36 2001
 * \brief  cdi_ipcress test function headers.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __cdi_ipcress_test_hh__
#define __cdi_ipcress_test_hh__

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <sstream>

namespace rtt_cdi_ipcress_test {

//---------------------------------------------------------------------------//
// COMPARISON FUNCTIONS USED IN IPCRESS OPACITY TESTS
//---------------------------------------------------------------------------//

template <typename temperatureType, typename densityType,
          typename testValueType, typename opacityClassType>
bool opacityAccessorPassed(rtt_dsxx::ScalarUnitTest &ut,
                           opacityClassType const spOpacity,
                           temperatureType const temperature,
                           densityType const density,
                           testValueType const tabulatedValue) {
  using std::ostringstream;

  // Interpolate the multigroup opacities.
  testValueType grayOpacity = spOpacity->getOpacity(temperature, density);

  // Make sure that the interpolated value matches previous
  // interpolations.

  if (rtt_dsxx::soft_equiv(grayOpacity, tabulatedValue)) {
    ostringstream message;
    message << spOpacity->getDataDescriptor()
            << " opacity computation was good for \n\t"
            << "\"" << spOpacity->getDataFilename() << "\" data.";
    ut.passes(message.str());
  } else {
    ostringstream message;
    message << spOpacity->getDataDescriptor()
            << " opacity value is out of spec. for \n\t"
            << "\"" << spOpacity->getDataFilename() << "\" data.";
    return ut.failure(message.str());
  }

  // If we get here then the test passed.
  return true;
}

//---------------------------------------------------------------------------//

template <typename opacityClassType>
void testTemperatureGridAccessor(rtt_dsxx::ScalarUnitTest &ut,
                                 opacityClassType const spOpacity) {
  using std::ostringstream;

  // Read the temperature grid from the IPCRESS file.
  std::vector<double> temps = spOpacity->getTemperatureGrid();

  // Verify that the size of the temperature grid looks right.  If
  // it is the right size then compare the temperature grid data to
  // the data specified when we created the IPCRESS file using TOPS.
  if (temps.size() == spOpacity->getNumTemperatures() && temps.size() == 3) {
    ostringstream message;
    message << "The number of temperature points found in the data\n\t"
            << "grid matches the number returned by the\n\t"
            << "getNumTemperatures() accessor.";
    ut.passes(message.str());

    // The grid specified by TOPS has 3 temperature points.
    std::vector<double> temps_ref(temps.size());
    temps_ref[0] = 0.1;
    temps_ref[1] = 1.0;
    temps_ref[2] = 10.0;

    // Compare the grids.
    if (rtt_dsxx::soft_equiv(temps, temps_ref))
      ut.passes("Temperature grid matches.");
    else
      ut.failure("Temperature grid did not match.");
  } else {
    ostringstream message;
    message << "The number of temperature points found in the data\n\t"
            << "grid does not match the number returned by the\n\t"
            << "getNumTemperatures() accessor. \n"
            << "Did not test the results returned by\n\t"
            << "getTemperatureGrid().";
    ut.failure(message.str());
  }
}

//---------------------------------------------------------------------------//

template <typename opacityClassType>
void testDensityGridAccessor(rtt_dsxx::ScalarUnitTest &ut,
                             opacityClassType const spOpacity) {
  using std::ostringstream;

  // Read the grid from the IPCRESS file.
  std::vector<double> density = spOpacity->getDensityGrid();

  // Verify that the size of the density grid looks right.  If
  // it is the right size then compare the density grid data to
  // the data specified when we created the IPCRESS file using TOPS.
  if (density.size() == 3 && density.size() == spOpacity->getNumDensities()) {
    ostringstream message;
    message << "The number of density points found in the data\n\t"
            << "grid matches the number returned by the\n\t"
            << "getNumDensities() accessor.";
    ut.passes(message.str());

    // The grid specified by TOPS has 3 density points
    std::vector<double> density_ref(density.size());
    density_ref[0] = 0.1;
    density_ref[1] = 0.5;
    density_ref[2] = 1.0;

    // Compare the grids.
    if (rtt_dsxx::soft_equiv(density, density_ref))
      ut.passes("Density grid matches.");
    else
      ut.failure("Density grid did not match.");
  } else {
    ostringstream message;
    message << "The number of density points found in the data\n\t"
            << "grid does not match the number returned by the\n\t"
            << "getNumDensities() accessor. \n"
            << "Did not test the results returned by\n\t"
            << "getDensityGrid().";
    ut.failure(message.str());
  }
}

//---------------------------------------------------------------------------//

template <typename opacityClassType>
void testEnergyBoundaryAccessor(rtt_dsxx::ScalarUnitTest &ut,
                                opacityClassType const spOpacity) {
  using std::ostringstream;

  // Read the grid from the IPCRESS file.
  std::vector<double> ebounds = spOpacity->getGroupBoundaries();

  // Verify that the size of the group boundary grid looks right.  If
  // it is the right size then compare the energy groups grid data to
  // the data specified when we created the IPCRESS file using TOPS.
  if (ebounds.size() == 13 &&
      ebounds.size() == spOpacity->getNumGroupBoundaries()) {
    ostringstream message;
    message << "The number of energy boundary points found in the data\n\t"
            << "grid matches the number returned by the\n\t"
            << "getNumGroupBoundaries() accessor.";
    ut.passes(message.str());

    // The grid specified by TOPS has 13 energy boundaries.
    std::vector<double> ebounds_ref(ebounds.size());
    ebounds_ref[0] = 0.01;
    ebounds_ref[1] = 0.03;
    ebounds_ref[2] = 0.07;
    ebounds_ref[3] = 0.1;
    ebounds_ref[4] = 0.3;
    ebounds_ref[5] = 0.7;
    ebounds_ref[6] = 1.0;
    ebounds_ref[7] = 3.0;
    ebounds_ref[8] = 7.0;
    ebounds_ref[9] = 10.0;
    ebounds_ref[10] = 30.0;
    ebounds_ref[11] = 70.0;
    ebounds_ref[12] = 100.0;

    // Compare the grids.
    if (rtt_dsxx::soft_equiv(ebounds, ebounds_ref))
      ut.passes("Energy group boundary grid matches.");
    else
      ut.failure("Energy group boundary grid did not match.");
  } else {
    ostringstream message;
    message << "The number of energy boundary points found in the data\n\t"
            << "grid does not match the number returned by the\n\t"
            << "get NumGroupBoundaries() accessor. \n"
            << "Did not test the results returned by\n\t"
            << "getGroupBoundaries().";
    ut.failure(message.str());
  }
}

} // end namespace rtt_cdi_ipcress_test

#endif // __cdi_ipcress_test_hh__

//---------------------------------------------------------------------------//
// end of cdi_ipcress/test/cdi_ipcress_test.hh
//---------------------------------------------------------------------------//
