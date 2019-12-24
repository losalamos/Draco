//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_CPEloss/test/tstTabular_CP_Eloss.cc
 * \author Ben R. Ryan
 * \date   2019 Nov 18
 * \brief  Tabular_CP_Eloss test
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "cdi_CPEloss/Tabular_CP_Eloss.hh"
#include "cdi/CDI.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/dbc.hh"
#include <sstream>

using rtt_cdi::CDI;
using rtt_cdi_cpeloss::Tabular_CP_Eloss;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void dedx_table_test(rtt_dsxx::UnitTest &ut) {

  // Datatable filename
  std::string filename_in = ut.getTestSourcePath() + "001-H-001";
  // Alpha particle target
  rtt_cdi::CParticle target_in(1002, 3.34447643e-24);
  // Proton projectile
  rtt_cdi::CParticle projectile_in(1001, 1.6726219e-24);

  Tabular_CP_Eloss eloss_mod(filename_in, target_in, projectile_in);

  // Model type better be tabular:
  FAIL_IF_NOT(eloss_mod.getModelType() == rtt_cdi::CPModelType::TABULAR_ETYPE);

  // Tabular data
  FAIL_IF_NOT(eloss_mod.data_in_tabular_form());

  // Check that grid accessors agree
  FAIL_IF_NOT(eloss_mod.getTemperatureGrid().size() ==
              eloss_mod.getNumTemperatures());
  FAIL_IF_NOT(eloss_mod.getDensityGrid().size() == eloss_mod.getNumDensities());
  FAIL_IF_NOT(eloss_mod.getEnergyGrid().size() == eloss_mod.getNumEnergies());

  // Check that class grid dimensions match file
  FAIL_IF_NOT(eloss_mod.getNumEnergies() == 4);
  FAIL_IF_NOT(eloss_mod.getNumDensities() == 5);
  FAIL_IF_NOT(eloss_mod.getNumTemperatures() == 6);

  // Data file should not be empty
  FAIL_IF(eloss_mod.getDataFilename().empty());

  // Throw error if trying to read table out of bounds
  try {
    eloss_mod.getEloss(0, 0, 0);
    FAILMSG("Require failed to fire with out-of-bounds table access.");
  } catch (rtt_dsxx::assertion &) {
    PASSMSG("Require fired for out-of-bounds table access.");
  }

  // Get eloss value for almost first (1,1,1) grid point
  {
    double energy = 1.384272;
    double density = 0.01314446475045816;
    double temperature = 1.9783097972766984e-05;
    printf("dedx: %28.18e\n", eloss_mod.getEloss(temperature, density, energy));
    printf("equiv? %i\n", (rtt_dsxx::soft_equiv(                                            
        eloss_mod.getEloss(temperature, density, energy), 2.652631e+03, 1.e-8)));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(
        //eloss_mod.getEloss(temperature, density, energy), 131.2710473, 1.e-8));
        eloss_mod.getEloss(temperature, density, energy), 2.652630801632265957e+03, 1.e-8));
  }

  // Get eloss value for almost last (2,3,4) grid point
  {
    double energy = 1.384273e+02;
    double density = 335853.3534222285;
    double temperature = 505.47663380532043;
    printf("dedx: %28.18e\n", eloss_mod.getEloss(temperature, density, energy));
    FAIL_IF_NOT(rtt_dsxx::soft_equiv(
        //eloss_mod.getEloss(temperature, density, energy), 40747641.884, 1.e-8));
        eloss_mod.getEloss(temperature, density, energy), 6.777724392862609375e+12, 1.e-8));
  }

  // Get eloss value for a point between grid points (45.5,13.5,13.5, i.e. requiring linear interpolation)
  // Get eloss value for a point between grid points (1.5,2.5,3.5, i.e. requiring linear interpolation)
  {
    double energy = 3.2272259644150236;
    double density = 66.4425508644999;
    double temperature = 0.09999946884616503;
    printf("dedx: %28.18e\n", eloss_mod.getEloss(temperature, density, energy));
    FAIL_IF_NOT(
        rtt_dsxx::soft_equiv(eloss_mod.getEloss(temperature, density, energy),
                             //155837.597258, 1.e-8));
                             3.125994000940805674e+07, 1.e-8));
  }

  if (ut.numFails == 0) {
    PASSMSG("Tabular_CP_Eloss test passes.");
  } else {
    FAILMSG("Tabular_CP_Eloss test fails.");
  }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    dedx_table_test(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tst_Tabular_CP_Eloss.cc
//---------------------------------------------------------------------------//
