//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_ipcress/test/tIpcressOpacity.cc
 * \author Thomas M. Evans
 * \date   Fri Oct 26 10:50:44 2001
 * \brief
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "cdi_ipcress_test.hh"
#include "cdi/OpacityCommon.hh"
#include "cdi_ipcress/IpcressFile.hh"
#include "cdi_ipcress/IpcressGrayOpacity.hh"
#include "cdi_ipcress/IpcressMultigroupOpacity.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <algorithm>

using namespace std;

using rtt_cdi_ipcress::IpcressGrayOpacity;
using rtt_cdi_ipcress::IpcressMultigroupOpacity;
using rtt_cdi_ipcress::IpcressFile;
using rtt_cdi::GrayOpacity;
using rtt_cdi::MultigroupOpacity;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void file_check_Al_BeCu(rtt_dsxx::ScalarUnitTest &ut) {
  // Ipcress data filename (IPCRESS format required)
  string op_data_file = "Al_BeCu.ipcress";

  // ------------------------- //
  // Create IpcressFile object //
  // ------------------------- //

  {
    shared_ptr<IpcressFile> spIF;

    // Attempt to instantiate the object.
    try {
      spIF.reset(new IpcressFile(op_data_file));
    } catch (rtt_dsxx::assertion const &error) {
      cout << "While testing tIpcressOpacity, " << error.what() << endl;
      return;
    }

    // If we make it here then spIF was successfully instantiated.
    PASSMSG(string("shared_ptr to new IpcressFile object created for ") +
            string("Al_BeCu.ipcress data."));

    // Test the IpcressFile object.
    if (spIF->getDataFilename() == op_data_file) {
      ostringstream message;
      message << "IpcressFile object is now linked to the "
              << "Al_BeCu.ipcress data file.";
      PASSMSG(message.str());
    } else {
      ostringstream message;
      message << "IpcressFile object failed to link itself to the "
              << "Al_BeCu.ipcress  data file.";
      FAILMSG(message.str());
    }

    if (spIF->getNumMaterials() == 2) {
      ostringstream message;
      message << "The correct number of materials was found in the "
              << "Al_BeCu.ipcress data file.";
      PASSMSG(message.str());
    } else {
      ostringstream message;
      message << "spIF did not find the correct number of materials "
              << "in the Al_BeCu.ipcress data file.";
      FAILMSG(message.str());
    }
  }

  // ---------------------- //
  // Create Opacity object. //
  // ---------------------- //

  // Material identifier.  This data file has two materials: Al and
  // BeCu.  Al has the id tag "10001".
  size_t const matid(10001);

  // Try to instantiate the Opacity object. (Rosseland, Gray Total
  // for material 10001 in the IPCRESS file pointed to by spIF).
  shared_ptr<GrayOpacity> spOp_Al_rgt;

  try {
    shared_ptr<IpcressFile> spIF(new IpcressFile(op_data_file));
    spOp_Al_rgt.reset(new IpcressGrayOpacity(spIF, matid, rtt_cdi::ROSSELAND,
                                             rtt_cdi::TOTAL));
    // spIF goes out of scope
  } catch (rtt_dsxx::assertion const &excpt) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for "
            << "Al_BeCu.ipcress data." << endl
            << "\t" << excpt.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  PASSMSG("shared_ptr to new Opacity object created for Al_BeCu.ipcress data.");

  // ----------------- //
  // Gray Opacity Test //
  // ----------------- //

  double temperature = 0.1;                         // keV
  double density = 27.0;                            // g/cm^3
  double tabulatedGrayOpacity = 4271.7041147070677; // cm^2/g

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(
          ut, spOp_Al_rgt, temperature, density, tabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // Check accessor functions

  rtt_cdi::OpacityModelType omt(spOp_Al_rgt->getOpacityModelType());
  if (omt == rtt_cdi::IPCRESS_TYPE)
    PASSMSG("OpacityModelType() returned expected value.");
  else
    FAILMSG("OpacityModelType() did not return the expected value.");

  string edp(spOp_Al_rgt->getEnergyPolicyDescriptor());
  if (edp == string("gray"))
    PASSMSG("EDP = gray");
  else
    FAILMSG("EDP != gray");

  if (!spOp_Al_rgt->data_in_tabular_form())
    ITFAILS;

  size_t nd(spOp_Al_rgt->getNumDensities());
  size_t nt(spOp_Al_rgt->getNumTemperatures());
  if (nd != 5)
    FAILMSG("Found wrong number of density values.");
  if (nt != 10)
    FAILMSG("Found wrong number of temperature values.");

  vector<double> densGrid(spOp_Al_rgt->getDensityGrid());
  vector<double> tempGrid(spOp_Al_rgt->getTemperatureGrid());
  if (densGrid.size() != nd)
    ITFAILS;
  if (tempGrid.size() != nt)
    ITFAILS;

  double expected_densGrid[] = {0.01, 0.1, 1.0, 10.0, 100.0};
  double expected_tempGrid[] = {0.0005, 0.0015, 0.004, 0.0125, 0.04,
                                0.125,  0.4,    1.25,  4,      15};

  for (size_t i = 0; i < densGrid.size(); ++i)
    if (!soft_equiv(densGrid[i], expected_densGrid[i]))
      ITFAILS;
  for (size_t i = 0; i < tempGrid.size(); ++i)
    if (!soft_equiv(tempGrid[i], expected_tempGrid[i]))
      ITFAILS;

  // --------------- //
  // MG Opacity test //
  // --------------- //

  // Create a Multigroup Rosseland Total Opacity object (again for Al).
  shared_ptr<MultigroupOpacity> spOp_Al_rtmg;

  // Try to instantiate the Opacity object.
  try {
    shared_ptr<IpcressFile> spIF(new IpcressFile(op_data_file));
    spOp_Al_rtmg.reset(new IpcressMultigroupOpacity(
        spIF, matid, rtt_cdi::ROSSELAND, rtt_cdi::TOTAL));
    // spIF goes out of scope
  } catch (rtt_dsxx::assertion const &excpt) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for "
            << "Al_BeCu.ipcress data." << endl
            << "\t" << excpt.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // Setup the test point.
  temperature = 0.01; // keV
  density = 2.0;      // g/cm^3

  // Check accessor functions

  omt = spOp_Al_rtmg->getOpacityModelType();
  if (omt == rtt_cdi::IPCRESS_TYPE)
    PASSMSG("OpacityModelType() returned expected value.");
  else
    FAILMSG("OpacityModelType() did not return the expected value.");

  edp = spOp_Al_rtmg->getEnergyPolicyDescriptor();
  if (edp == string("mg"))
    PASSMSG("EDP = mg");
  else
    FAILMSG("EDP != mg");

  if (!spOp_Al_rtmg->data_in_tabular_form())
    ITFAILS;

  size_t numGroups = 33;
  if (spOp_Al_rtmg->getNumGroups() != numGroups)
    ITFAILS;

  // The solution to compare against:
  double tabulatedMGOpacityArray[] = {
      2.4935245299837247e+08, 2.6666789027326573e+04,
      1.6270621515227660e+04, 1.7634711671468522e+04,
      4.4999455359684442e+04, 9.9917674335613032e+04,
      8.3261383385464113e+04, 5.9742975304886764e+04,
      4.0373209722602740e+04, 2.6156503146710609e+04,
      1.6356701105166874e+04, 1.0007184686170869e+04,
      5.9763667878215247e+03, 3.5203912630050986e+03,
      2.0765528559140448e+03, 6.8550529299142445e+03,
      4.1257095227186965e+03, 2.4199006949490426e+03,
      1.3894677080938793e+03, 7.9046985091966621e+02,
      4.4088463936537232e+02, 2.4514360684176387e+02,
      1.3541656611912146e+02, 7.1828886317050177e+01,
      3.9793827527329107e+01, 2.3312673181867030e+01,
      1.4879458895157605e+01, 1.0862672679200283e+01,
      9.0590676798691288e+00, 8.2841367649864175e+00,
      7.3809286930540363e+00, 7.1057875403123791e+00,
      6.8907716134926735e+00}; // KeV, numGroups entries.

  vector<double> tabulatedMGOpacity(numGroups);
  copy(tabulatedMGOpacityArray, tabulatedMGOpacityArray + numGroups,
       tabulatedMGOpacity.begin());

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(
          ut, spOp_Al_rtmg, temperature, density, tabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }
}

//---------------------------------------------------------------------------//

void file_check_analytic(rtt_dsxx::ScalarUnitTest &ut) {
  // ----------------------------------------------------------------
  // The Opacities in this file are computed from the following
  // analytic formula:
  //     opacity = rho * T^4,
  // rho is the density and T is the temperature.
  //
  // The grid in this data file has the following structure:
  //    T   = { 0.1, 1.0, 10.0 } keV.
  //    rho = { 0.1, 0.5, 1.0 } g/cm^3
  //    E_bounds = { 0.01, 0.03, 0.07, 0.1, 0.3, 0.7, 1.0, 3.0, 7.0
  //                 10.0, 30.0, 70.0 100.0 } keV.
  //-----------------------------------------------------------------

  // Ipcress data filename (IPCRESS format required)
  string op_data_file = "analyticOpacities.ipcress";

  // ------------------------- //
  // Create IpcressFile object //
  // ------------------------- //

  // Create a smart pointer to a IpcressFile object
  shared_ptr<IpcressFile> spGFAnalytic;

  // Try to instantiate the object.
  try {
    spGFAnalytic.reset(new rtt_cdi_ipcress::IpcressFile(op_data_file));
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    FAILMSG(error.what());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we make it here then spGFAnalytic was successfully instantiated.
  PASSMSG("shared_ptr to new IpcressFile object created (spGFAnalytic).");

  // Test the IpcressFile object.
  if (spGFAnalytic->getDataFilename() == op_data_file) {
    ostringstream message;
    message << "IpcressFile object is now linked to the data file.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "IpcressFile object failed to link itself "
            << "to the data file.";
    FAILMSG(message.str());
  }

  if (spGFAnalytic->getNumMaterials() == 1) {
    ostringstream message;
    message << "The correct number of materials was found "
            << "in the data file.";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << "spGFAnalytic did not find the correct number "
            << "of materials in the data file.";
    FAILMSG(message.str());
  }

  // --------------------- //
  // Create Opacity object //
  // --------------------- //

  // Create a smart pointer to an Opacity object.
  shared_ptr<GrayOpacity> spOp_Analytic_ragray;

  // material ID
  const int matid = 10001;

  // Try to instantiate the Opacity object.
  try {
    spOp_Analytic_ragray.reset(new IpcressGrayOpacity(
        spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION));
  } catch (rtt_dsxx::assertion const &error)
  // Alternatively, we could use:
  // catch ( rtt_cdi_ipcress::gkeysException GandError )
  // catch ( rtt_cdi_ipcress::gchgridsException GandError )
  // catch ( rtt_cdi_ipcress::ggetmgException GandError )
  // catch ( rtt_cdi_ipcress::ggetgrayException GandError )
  {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for "
            << "analyticOpacities.ipcress data." << endl
            << "\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  PASSMSG(string("shared_ptr to new Opacity object created for ") +
          string("analyticOpacities.ipcress."));

  // ----------------- //
  // Gray Opacity Test //
  // ----------------- //

  double temperature = 10.0; // keV
  double density = 1.0;      // g/cm^3
  double tabulatedGrayOpacity = density * pow(temperature, 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_ragray,
                                                   temperature, density,
                                                   tabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  //---------------- //
  // MG Opacity test //
  //---------------- //

  // Create a smart pointer to an Opacity object.
  shared_ptr<MultigroupOpacity> spOp_Analytic_ramg;

  // Try to instantiate the Opacity object.
  try {
    spOp_Analytic_ramg.reset(new IpcressMultigroupOpacity(
        spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION));
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for "
            << "analyticOpacities.ipcress data." << endl
            << "\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  PASSMSG(string("shared_ptr to new Opacity object created for ") +
          string("analyticOpacities.ipcress."));

  // Set up the new test problem.

  temperature = 0.3; // keV
  density = 0.7;     // g/cm^3

  // This is the solution we compare against.
  int numGroups = 12;
  vector<double> tabulatedMGOpacity(numGroups);
  for (int i = 0; i < numGroups; ++i)
    tabulatedMGOpacity[i] = density * pow(temperature, 4); // cm^2/gm

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(
          ut, spOp_Analytic_ramg, temperature, density, tabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------------------------------------------------ //
  // Test the Plank routines using analyticOpacities.ipcress data //
  // ------------------------------------------------------------ //

  // The Opacities in this file are computed from the following
  // analytic formula:
  //     opacity = rho * T^4,
  // rho is the density and T is the temperature.

  // spGFAnalytic already points to the correct file so we don't repeat the
  // coding.

  // Dito for spOpacityAnalytic.

  // ----------------- //
  // Gray Opacity Test //
  // ----------------- //

  // Create a smart pointer to an Opacity object.
  shared_ptr<GrayOpacity> spOp_Analytic_pgray;

  // Try to instantiate the Opacity object.
  try {
    spOp_Analytic_pgray.reset(new IpcressGrayOpacity(
        spGFAnalytic, matid, rtt_cdi::PLANCK, rtt_cdi::ABSORPTION));
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for "
            << "analyticOpacities.ipcress data." << endl
            << "\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  {
    ostringstream message;
    message << "shared_ptr to new Gray Plank Total Opacity object "
            << "created for analyticOpacities.ipcress.";
    PASSMSG(message.str());
  }

  // Setup the test problem.

  temperature = 3.0;                                     // keV
  density = 0.7;                                         // g/cm^3
  double tabulatedValue = density * pow(temperature, 4); // cm^2/g

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(
          ut, spOp_Analytic_pgray, temperature, density, tabulatedValue)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // --------------- //
  // MG Opacity test //
  // --------------- //

  // Create a smart pointer to an Opacity object.
  shared_ptr<MultigroupOpacity> spOp_Analytic_pmg;

  // Try to instantiate the Opacity object.
  try {
    spOp_Analytic_pmg.reset(new IpcressMultigroupOpacity(
        spGFAnalytic, matid, rtt_cdi::PLANCK, rtt_cdi::ABSORPTION));
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressOpacity object for "
            << "analyticOpacities.ipcress data." << endl
            << "\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  {
    ostringstream message;
    message << "shared_ptr to new Multigroup Plank Total Opacity object "
            << "created \n\t for \"analyticOpacities.ipcress.\"";
    PASSMSG(message.str());
  }

  // Setup the test problem.

  size_t ng = 12;
  tabulatedMGOpacity.resize(ng);
  temperature = 0.4; // keV
  density = 0.22;    // g/cm^3
  for (size_t ig = 0; ig < ng; ++ig)
    tabulatedMGOpacity[ig] = density * pow(temperature, 4); // cm^2/g

  // If this test fails then stop testing.
  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(
          ut, spOp_Analytic_pmg, temperature, density, tabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------------ //
  // Access temperature grid. //
  // ------------------------ //

  rtt_cdi_ipcress_test::testTemperatureGridAccessor(ut, spOp_Analytic_pmg);

  // ------------------------ //
  // Access the density grid. //
  // ------------------------ //

  rtt_cdi_ipcress_test::testDensityGridAccessor(ut, spOp_Analytic_pmg);

  // ----------------------------- //
  // Access the energy boundaries. //
  // ----------------------------- //

  rtt_cdi_ipcress_test::testEnergyBoundaryAccessor(ut, spOp_Analytic_pmg);

  // ------------------------------------------------------------ //
  // Test alternate (vector-based) accessors for getGrayRosseland //
  // ------------------------------------------------------------ //

  // ---------------------- //
  // Vector of temperatures //
  // ---------------------- //

  vector<double> vtemperature(2);
  vtemperature[0] = 0.5; // keV
  vtemperature[1] = 0.7; // keV
  density = 0.35;        // g/cm^3

  vector<double> vtabulatedGrayOpacity(vtemperature.size());
  for (size_t i = 0; i < vtabulatedGrayOpacity.size(); ++i)
    vtabulatedGrayOpacity[i] = density * pow(vtemperature[i], 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_ragray,
                                                   vtemperature, density,
                                                   vtabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ---------------------- //
  // Vector of densities    //
  // ---------------------- //

  vector<double> vdensity(3);
  temperature = 0.3; //keV
  vdensity[0] = 0.2; // g/cm^3
  vdensity[1] = 0.4; // g/cm^3
  vdensity[2] = 0.6; // g/cm^3

  vtabulatedGrayOpacity.resize(vdensity.size());
  for (size_t i = 0; i < vtabulatedGrayOpacity.size(); ++i)
    vtabulatedGrayOpacity[i] = vdensity[i] * pow(temperature, 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_ragray,
                                                   temperature, vdensity,
                                                   vtabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // -------------------------------------------------------- //
  // Test alternate (vector-based) accessors for getGrayPlank //
  // -------------------------------------------------------- //

  // ---------------------- //
  // Vector of temperatures //
  // ---------------------- //

  vtemperature.resize(2);
  vtemperature[0] = 0.5; // keV
  vtemperature[1] = 0.7; // keV
  density = 0.35;        // g/cm^3q

  vtabulatedGrayOpacity.resize(vtemperature.size());
  for (size_t i = 0; i < vtabulatedGrayOpacity.size(); ++i)
    vtabulatedGrayOpacity[i] = density * pow(vtemperature[i], 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_pgray,
                                                   vtemperature, density,
                                                   vtabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------- //
  // Vector of densities //
  // ------------------- //

  vdensity.resize(3);
  temperature = 0.3; //keV
  vdensity[0] = 0.2; // g/cm^3
  vdensity[1] = 0.4; // g/cm^3
  vdensity[2] = 0.6; // g/cm^3

  vtabulatedGrayOpacity.resize(vdensity.size());
  for (size_t i = 0; i < vtabulatedGrayOpacity.size(); ++i)
    vtabulatedGrayOpacity[i] = vdensity[i] * pow(temperature, 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(ut, spOp_Analytic_pgray,
                                                   temperature, vdensity,
                                                   vtabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ---------------------------------------------------------- //
  // Test alternate (vector-based) accessors for getMGRosseland //
  // ---------------------------------------------------------- //

  // ---------------------- //
  // Vector of temperatures //
  // ---------------------- //

  vtemperature.resize(2);
  vtemperature[0] = 0.5; // keV
  vtemperature[1] = 0.7; // keV
  density = 0.35;        // g/cm^3
  ng = spOp_Analytic_ramg->getNumGroupBoundaries() - 1;

  vector<vector<double>> vtabulatedMGOpacity(vtemperature.size());
  for (size_t i = 0; i < vtemperature.size(); ++i) {
    vtabulatedMGOpacity[i].resize(ng);
    for (size_t ig = 0; ig < ng; ++ig)
      vtabulatedMGOpacity[i][ig] = density * pow(vtemperature[i], 4);
  }

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(
          ut, spOp_Analytic_ramg, vtemperature, density, vtabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------- //
  // Vector of densities //
  // ------------------- //

  vdensity.resize(2);
  vdensity[0] = 0.3; // g/cm^3
  vdensity[1] = 0.7; // g/cm^3
  temperature = 7.0; // keV
  ng = spOp_Analytic_ramg->getNumGroupBoundaries() - 1;

  vtabulatedMGOpacity.resize(vdensity.size());
  for (size_t i = 0; i < vdensity.size(); ++i) {
    vtabulatedMGOpacity[i].resize(ng);
    for (size_t ig = 0; ig < ng; ++ig)
      vtabulatedMGOpacity[i][ig] = vdensity[i] * pow(temperature, 4);
  }

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(
          ut, spOp_Analytic_ramg, temperature, vdensity, vtabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------------------------------------------ //
  // Test alternate (vector-based) accessors for getMGPlank //
  // ------------------------------------------------------ //

  // ---------------------- //
  // Vector of temperatures //
  // ---------------------- //

  vtemperature.resize(2);
  vtemperature[0] = 0.5; // keV
  vtemperature[1] = 0.7; // keV
  density = 0.35;        // g/cm^3
  ng = spOp_Analytic_pmg->getNumGroupBoundaries() - 1;

  vtabulatedMGOpacity.resize(vtemperature.size());
  for (size_t i = 0; i < vtemperature.size(); ++i) {
    vtabulatedMGOpacity[i].resize(ng);
    for (size_t ig = 0; ig < ng; ++ig)
      vtabulatedMGOpacity[i][ig] = density * pow(vtemperature[i], 4);
  }

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(
          ut, spOp_Analytic_pmg, vtemperature, density, vtabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------- //
  // Vector of densities //
  // ------------------- //

  vdensity.resize(2);
  vdensity[0] = 0.3; // g/cm^3
  vdensity[1] = 0.7; // g/cm^3
  temperature = 7.0; // keV
  ng = spOp_Analytic_pmg->getNumGroupBoundaries() - 1;

  vtabulatedMGOpacity.resize(vdensity.size());
  for (size_t i = 0; i < vdensity.size(); ++i) {
    vtabulatedMGOpacity[i].resize(ng);
    for (size_t ig = 0; ig < ng; ++ig)
      vtabulatedMGOpacity[i][ig] = vdensity[i] * pow(temperature, 4);
  }

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(
          ut, spOp_Analytic_pmg, temperature, vdensity, vtabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }
}

//---------------------------------------------------------------------------//

void check_ipcress_stl_accessors(rtt_dsxx::ScalarUnitTest &ut) {
  // Ipcress data filename (IPCRESS format required)
  string op_data_file = "analyticOpacities.ipcress";

  // ------------------------- //
  // Create IpcressFile object //
  // ------------------------- //

  // Create a smart pointer to a IpcressFile object
  shared_ptr<IpcressFile> spGFAnalytic;

  // Try to instantiate the object.
  try {
    spGFAnalytic.reset(new rtt_cdi_ipcress::IpcressFile(op_data_file));
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    FAILMSG(error.what());
    FAILMSG("Aborting tests.");
    return;
  }

  // material ID
  const int matid = 10001;

  // -------------------------------------- //
  // Test the STL-like getOpacity accessor  //
  // Using const iterators for Gray objects //
  // -------------------------------------- //

  // These accessors are only available in IpcressOpacity objects so the
  // shared_ptr must be templated on IpcressGrayOpacity and not on
  // cdi/GrayOpacity.

  // Create a new smart pointer to a IpcressGrayOpacity object.
  shared_ptr<IpcressGrayOpacity> spGGOp_Analytic_ra;

  // try to instantiate the Opacity object.
  try {
    spGGOp_Analytic_ra.reset(new IpcressGrayOpacity(
        spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION));
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressGrayOpacity object "
            << "fpr \n\t analyticOpacityies.ipcress data (shared_ptr not "
            << "templated on cdi/GrayOpacity).\n\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  vector<double> vdensity;
  vector<double> vtemperature;
  vector<double> vtabulatedGrayOpacity;

  double density = 0.35;    // g/cm^3
  double temperature = 7.0; // kev

  // Setup the temperature and density parameters for this test.
  vdensity.resize(6);
  vtemperature.resize(6);

  // (temperature,density) tuples.

  vtemperature[0] = 0.5; // keV
  vdensity[0] = 0.2;     // g/cm^3

  vtemperature[1] = 0.7; // keV
  vdensity[1] = 0.2;     // g/cm^3

  vtemperature[2] = 0.5; // keV
  vdensity[2] = 0.4;     // g/cm^3

  vtemperature[3] = 0.7; // keV
  vdensity[3] = 0.4;     // g/cm^3

  vtemperature[4] = 0.5; // keV
  vdensity[4] = 0.6;     // g/cm^3

  vtemperature[5] = 0.7; // keV
  vdensity[5] = 0.6;     // g/cm^3

  // we want to test the const_iterator version of getOpacity() so
  // we need to create const vectors with the tuple data.
  const vector<double> cvdensity = vdensity;
  const vector<double> cvtemperature = vtemperature;

  int nt = cvtemperature.size();
  int nd = cvdensity.size();

  // Here is the reference solution
  vtabulatedGrayOpacity.resize(nt);
  for (int i = 0; i < nt; ++i)
    vtabulatedGrayOpacity[i] = cvdensity[i] * pow(cvtemperature[i], 4);

  // Here is the solution from Ipcress
  vector<double> graOp(nt);
  spGGOp_Analytic_ra->getOpacity(cvtemperature.begin(), cvtemperature.end(),
                                 cvdensity.begin(), cvdensity.end(),
                                 graOp.begin());
  if (rtt_dsxx::soft_equiv(graOp, vtabulatedGrayOpacity)) {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor()
            << " opacity computation was good for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (const-iterator accessor, temp x density).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor()
            << " opacity value is out of spec. for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    FAILMSG(message.str());
  }

  // ------------------------------------- //
  // Test the STL-like getOpacity accessor //
  // Using non-const iterator              //
  // ------------------------------------- //

  spGGOp_Analytic_ra->getOpacity(vtemperature.begin(), vtemperature.end(),
                                 vdensity.begin(), vdensity.end(),
                                 graOp.begin());
  if (rtt_dsxx::soft_equiv(graOp, vtabulatedGrayOpacity)) {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor()
            << " opacity computation was good for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor()
            << " opacity value is out of spec. for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    FAILMSG(message.str());
  }

  // ------------------------------------- //
  // Test the STL-like getOpacity accessor //
  // const iterator (temperature only)     //
  // ------------------------------------- //

  graOp.resize(nt);
  vtabulatedGrayOpacity.resize(nt);
  for (int it = 0; it < nt; ++it)
    vtabulatedGrayOpacity[it] = density * pow(vtemperature[it], 4);

  spGGOp_Analytic_ra->getOpacity(cvtemperature.begin(), cvtemperature.end(),
                                 density, graOp.begin());
  if (rtt_dsxx::soft_equiv(graOp, vtabulatedGrayOpacity)) {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor()
            << " opacity computation was good for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (const iterator accessor, vtemps).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor()
            << " opacity value is out of spec. for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (const iterator accessor, vtemps).";
    FAILMSG(message.str());
  }

  // ------------------------------------- //
  // Test the STL-like getOpacity accessor //
  // const iterator ( density only)        //
  // ------------------------------------- //

  graOp.resize(nd);
  vtabulatedGrayOpacity.resize(nd);
  for (int id = 0; id < nd; ++id)
    vtabulatedGrayOpacity[id] = vdensity[id] * pow(temperature, 4);

  spGGOp_Analytic_ra->getOpacity(temperature, cvdensity.begin(),
                                 cvdensity.end(), graOp.begin());
  if (rtt_dsxx::soft_equiv(graOp, vtabulatedGrayOpacity)) {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor()
            << " opacity computation was good for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (const iterator accessor, vdensity).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGGOp_Analytic_ra->getDataDescriptor()
            << " opacity value is out of spec. for \n\t"
            << spGGOp_Analytic_ra->getDataFilename()
            << " (const iterator accessor, vdensity).";
    FAILMSG(message.str());
  }

  // -------------------------------------- //
  // Test the STL-like getOpacity accessor  //
  // Using const iterators for MG objects   //
  // -------------------------------------- //

  // These accessors are only available in IpcressOpacity objects
  // so the shared_ptr must be templated on IpcressMultigroupOpacity and not on
  // cdi/MultigroupOpacity.

  // Create a new smart pointer to a IpcressGrayOpacity object.
  shared_ptr<IpcressMultigroupOpacity> spGMGOp_Analytic_ra;

  // try to instantiate the Opacity object.
  try {
    spGMGOp_Analytic_ra.reset(new IpcressMultigroupOpacity(
        spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION));
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to new IpcressGrayOpacity "
            << "object for \n\t analyticOpacities.ipcress data "
            << "(shared_ptr not templated on cdi/GrayOpacity)." << endl
            << "\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // If we get here then the object was successfully instantiated.
  PASSMSG(string("shared_ptr to new Opacity object created for ") +
          string("analyticOpacities.ipcress."));

  // Here is the reference solution
  int ng = spGMGOp_Analytic_ra->getNumGroupBoundaries() - 1;
  vector<double> vtabulatedOpacity(ng * nt);

  for (int i = 0; i < nt; ++i)
    for (int ig = 0; ig < ng; ++ig)
      vtabulatedOpacity[i * ng + ig] = cvdensity[i] * pow(cvtemperature[i], 4);

  // Here is the solution from Ipcress
  vector<double> mgOp(nt * ng);
  spGMGOp_Analytic_ra->getOpacity(cvtemperature.begin(), cvtemperature.end(),
                                  cvdensity.begin(), cvdensity.end(),
                                  mgOp.begin());

  if (rtt_dsxx::soft_equiv(mgOp, vtabulatedOpacity)) {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor()
            << " opacity computation was good for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (const-iterator accessor, temp x density).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor()
            << " opacity value is out of spec. for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    FAILMSG(message.str());
  }

  // --------------------------------------- //
  // Test the STL-like getOpacity accessor   //
  // Using non-const iterator for MG objects //
  // --------------------------------------- //

  // clear old data
  for (int i = 0; i < nt * ng; ++i)
    mgOp[i] = 0.0;

  // use Ipcress to obtain new data
  spGMGOp_Analytic_ra->getOpacity(vtemperature.begin(), vtemperature.end(),
                                  vdensity.begin(), vdensity.end(),
                                  mgOp.begin());

  // compare the results to the reference solution and report our findings.
  if (rtt_dsxx::soft_equiv(mgOp, vtabulatedOpacity)) {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor()
            << " opacity computation was good for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor()
            << " opacity value is out of spec. for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (non-const-iterator accessor, temp x density).";
    FAILMSG(message.str());
  }

  // ------------------------------------------------ //
  // Test the STL-like getOpacity accessor            //
  // const iterator (temperature only) for MG data    //
  // ------------------------------------------------ //

  // clear old data
  for (int i = 0; i < nt * ng; ++i)
    mgOp[i] = 0.0;

  // Calculate the reference solution.
  for (int it = 0; it < nt; ++it)
    for (int ig = 0; ig < ng; ++ig)
      vtabulatedOpacity[it * ng + ig] = density * pow(vtemperature[it], 4);

  // Obtain new solution
  spGMGOp_Analytic_ra->getOpacity(cvtemperature.begin(), cvtemperature.end(),
                                  density, mgOp.begin());

  // Compare solutions and report the results.
  if (rtt_dsxx::soft_equiv(mgOp, vtabulatedOpacity)) {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor()
            << " opacity computation was good for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (const iterator accessor, vtemps).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor()
            << " opacity value is out of spec. for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (const iterator accessor, vtemps).";
    FAILMSG(message.str());
  }

  // ------------------------------------------ //
  // Test the STL-like getOpacity accessor      //
  // const iterator ( density only) for MG data //
  // ------------------------------------------ //

  // clear old data
  for (int i = 0; i < nd * ng; ++i)
    mgOp[i] = 0.0;

  // Calculate the reference solution.
  for (int id = 0; id < nd; ++id)
    for (int ig = 0; ig < ng; ++ig)
      vtabulatedOpacity[id * ng + ig] = vdensity[id] * pow(temperature, 4);

  // Obtain new solution
  spGMGOp_Analytic_ra->getOpacity(temperature, cvdensity.begin(),
                                  cvdensity.end(), mgOp.begin());

  // Compare solutions and report the results.
  if (rtt_dsxx::soft_equiv(mgOp, vtabulatedOpacity)) {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor()
            << " opacity computation was good for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (const iterator accessor, vdensity).";
    PASSMSG(message.str());
  } else {
    ostringstream message;
    message << spGMGOp_Analytic_ra->getDataDescriptor()
            << " opacity value is out of spec. for \n\t"
            << spGMGOp_Analytic_ra->getDataFilename()
            << " (const iterator accessor, vdensity).";
    FAILMSG(message.str());
  }
}

//---------------------------------------------------------------------------//

void gray_opacity_packing_test(rtt_dsxx::ScalarUnitTest &ut) {
  vector<char> packed;

  {
    // Ipcress data filename (IPCRESS format required)
    string op_data_file = "analyticOpacities.ipcress";

    // ------------------------- //
    // Create IpcressFile object //
    // ------------------------- //

    // Create a smart pointer to a IpcressFile object
    shared_ptr<IpcressFile> spGFAnalytic;

    // Try to instantiate the object.
    try {
      spGFAnalytic.reset(new rtt_cdi_ipcress::IpcressFile(op_data_file));
    } catch (rtt_dsxx::assertion const &error) {
      ostringstream message;
      FAILMSG(error.what());
      FAILMSG("Aborting tests.");
      return;
    }

    // Create a smart pointer to an Opacity object.
    shared_ptr<GrayOpacity> spOp_Analytic_ragray;

    // material ID
    int const matid = 10001;

    spOp_Analytic_ragray.reset(new IpcressGrayOpacity(
        spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION));

    // pack up the opacity
    packed = spOp_Analytic_ragray->pack();
  }

  // make a new IpcressGrayOpacity from packed data
  shared_ptr<GrayOpacity> unpacked_opacity;

  // Try to instantiate the Opacity object.
  try {
    unpacked_opacity.reset(new IpcressGrayOpacity(packed));
  } catch (rtt_dsxx::assertion const &error) {
    ostringstream message;
    message << "Failed to create shared_ptr to unpacked IpcressOpacity object "
            << "for analyticOpacities.ipcress data." << endl
            << "\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // some simple tests
  if (unpacked_opacity->getDataFilename() != "analyticOpacities.ipcress")
    ITFAILS;

  if (unpacked_opacity->getReactionType() != rtt_cdi::ABSORPTION)
    ITFAILS;
  if (unpacked_opacity->getModelType() != rtt_cdi::ROSSELAND)
    ITFAILS;

  // ----------------- //
  // Gray Opacity Test //
  // ----------------- //

  double temperature = 10.0; // keV
  double density = 1.0;      // g/cm^3
  double tabulatedGrayOpacity = density * pow(temperature, 4);

  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(
          ut, unpacked_opacity, temperature, density, tabulatedGrayOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // try to unpack gray opacity as multigroup opacity
  bool caught = false;
  try {
    shared_ptr<MultigroupOpacity> opacity(new IpcressMultigroupOpacity(packed));
  } catch (rtt_dsxx::assertion const &error) {
    caught = true;
    ostringstream message;
    message << "Good, we caught the following assertion, \n" << error.what();
    PASSMSG(message.str());
  }
  if (!caught) {
    FAILMSG("Failed to catch an illegal packing asserion.");
  }
}

//---------------------------------------------------------------------------//

void mg_opacity_packing_test(rtt_dsxx::ScalarUnitTest &ut) {
  vector<char> packed;

  {
    // Ipcress data filename (IPCRESS format required)
    string op_data_file = "analyticOpacities.ipcress";

    // ------------------------- //
    // Create IpcressFile object //
    // ------------------------- //

    // Create a smart pointer to a IpcressFile object
    shared_ptr<IpcressFile> spGFAnalytic;

    // Try to instantiate the object.
    try {
      spGFAnalytic.reset(new rtt_cdi_ipcress::IpcressFile(op_data_file));
    } catch (rtt_dsxx::assertion const &error) {
      ostringstream message;
      FAILMSG(error.what());
      FAILMSG("Aborting tests.");
      return;
    }

    // material ID
    const int matid = 10001;

    //---------------- //
    // MG Opacity test //
    //---------------- //

    // Create a smart pointer to an Opacity object.
    shared_ptr<MultigroupOpacity> spOp_Analytic_pmg;

    spOp_Analytic_pmg.reset(new IpcressMultigroupOpacity(
        spGFAnalytic, matid, rtt_cdi::PLANCK, rtt_cdi::ABSORPTION));

    packed = spOp_Analytic_pmg->pack();
  }

  // make a new IpcressGrayOpacity from packed data
  shared_ptr<MultigroupOpacity> unpacked_opacity;

  // Try to instantiate the Opacity object.
  try {
    unpacked_opacity.reset(new IpcressMultigroupOpacity(packed));
  } catch (rtt_dsxx::assertion const &error)
  // Alternatively, we could use:
  // catch ( rtt_cdi_ipcress::gkeysException GandError )
  // catch ( rtt_cdi_ipcress::gchgridsException GandError )
  // catch ( rtt_cdi_ipcress::ggetmgException GandError )
  // catch ( rtt_cdi_ipcress::ggetgrayException GandError )
  {
    ostringstream message;
    message << "Failed to create shared_ptr to unpacked IpcressOpacity object "
            << "for analyticOpacities.ipcress data." << endl
            << "\t" << error.what();
    FAILMSG(message.str());
    FAILMSG("Aborting tests.");
    return;
  }

  // some simple tests
  if (unpacked_opacity->getDataFilename() != "analyticOpacities.ipcress")
    ITFAILS;

  if (unpacked_opacity->getReactionType() != rtt_cdi::ABSORPTION)
    ITFAILS;
  if (unpacked_opacity->getModelType() != rtt_cdi::PLANCK)
    ITFAILS;

  // Setup the test problem.

  int ng = 12;
  vector<double> tabulatedMGOpacity(ng);
  double temperature = 0.4; // keV
  double density = 0.22;    // g/cm^3
  for (int ig = 0; ig < ng; ++ig)
    tabulatedMGOpacity[ig] = density * pow(temperature, 4); // cm^2/g

  // If this test fails then stop testing.
  if (!rtt_cdi_ipcress_test::opacityAccessorPassed(
          ut, unpacked_opacity, temperature, density, tabulatedMGOpacity)) {
    FAILMSG("Aborting tests.");
    return;
  }

  // ------------------------ //
  // Access temperature grid. //
  // ------------------------ //

  rtt_cdi_ipcress_test::testTemperatureGridAccessor(ut, unpacked_opacity);

  // ------------------------ //
  // Access the density grid. //
  // ------------------------ //

  rtt_cdi_ipcress_test::testDensityGridAccessor(ut, unpacked_opacity);

  // ----------------------------- //
  // Access the energy boundaries. //
  // ----------------------------- //

  rtt_cdi_ipcress_test::testEnergyBoundaryAccessor(ut, unpacked_opacity);

  // try to unpack multigroup as gray opacity
  bool caught = false;
  try {
    shared_ptr<GrayOpacity> opacity(new IpcressGrayOpacity(packed));
  } catch (rtt_dsxx::assertion const &error) {
    caught = true;
    ostringstream message;
    message << "Good, we caught the following assertion, \n" << error.what();
    PASSMSG(message.str());
  }
  if (!caught) {
    FAILMSG("Failed to catch an illegal packing asserion.");
  }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    file_check_Al_BeCu(ut);
    file_check_analytic(ut);
    check_ipcress_stl_accessors(ut);

    gray_opacity_packing_test(ut);
    mg_opacity_packing_test(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tIpcressOpacity.cc
//---------------------------------------------------------------------------//
