//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/test/tstAnalytic_CP_Eloss.cc
 * \author Kendra P. Long
 * \date   Fri Aug  2 14:28:08 2019
 * \brief  Analytic_CP_Eloss test.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "cdi_analytic_test.hh"
#include "cdi/CDI.hh"
#include "cdi_analytic/Analytic_CP_Eloss.hh"
#include "cdi_analytic/Analytic_Models.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/dbc.hh"
#include <sstream>

using namespace std;

using rtt_cdi::CDI;
using rtt_cdi_analytic::Analytic_CP_Eloss;
using rtt_cdi_analytic::Analytic_Eloss_Model;
//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void KP_alpha_test(rtt_dsxx::UnitTest &ut) {

  // Deuterium:
  int32_t target_in = 1002;
  // Alpha particle:
  int32_t projectile_in = 2004;

  std::shared_ptr<Analytic_Eloss_Model> model_in(
      new rtt_cdi_analytic::Analytic_KP_Alpha_Eloss_Model());

  Analytic_CP_Eloss eloss_mod(model_in, target_in, projectile_in);

  // Check that basic accessors return correct result:
  // Model type better be analytic:
  FAIL_IF_NOT(eloss_mod.getModelType() == rtt_cdi::CPModelType::ANALYTIC_ETYPE);
  FAIL_IF_NOT(eloss_mod.getModel() == rtt_cdi::CPModel::ELOSS);

  // NOT tabular data
  FAIL_IF(eloss_mod.data_in_tabular_form());

  // All sizes should be 0 (again, not tabular data)
  FAIL_IF_NOT(eloss_mod.getTemperatureGrid().size() ==
              eloss_mod.getNumTemperatures());
  FAIL_IF_NOT(eloss_mod.getDensityGrid().size() == eloss_mod.getNumDensities());
  FAIL_IF_NOT(eloss_mod.getEnergyGrid().size() == eloss_mod.getNumEnergies());

  // Data file name should be an empty string:
  FAIL_IF_NOT(eloss_mod.getDataFilename().empty());

  // Check that accessors return the correct target and projectile:
  FAIL_IF_NOT(target_in == eloss_mod.getTargetZAID());
  FAIL_IF_NOT(projectile_in == eloss_mod.getProjectileZAID());

  // Get eloss values for some sample data:
  {
    double T = 1.0;    // keV
    double rho = 10.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 212.287, 1.0e-3));
  }

  // Point near the maximum on rho in [0, 20] g/cc, T in [0, 100] keV:
  {
    double T = 15.0;   // keV
    double rho = 20.0; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 78.5845, 1.0e-4));
  }

  // Point where the range fit goes negative -- eloss should return 1.0e25:
  // (equivalent to the particle range limiting to a very small number)
  {
    double T = 90.0;   // keV
    double rho = 17.5; // g / cc
    double vel0 = 1.0; // cm / shk

    double eloss_coeff = eloss_mod.getEloss(T, rho, vel0);

    FAIL_IF_NOT(rtt_dsxx::soft_equiv(eloss_coeff, 1.0e25));
  }

  if (ut.numFails == 0)
    PASSMSG("Analytic_CP_Eloss test passes.");
  else
    FAILMSG("Analytic_CP_Eloss test fails.");

  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    KP_alpha_test(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstAnalytic_CP_Eloss.cc
//---------------------------------------------------------------------------//
