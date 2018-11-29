//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/test/tstAnalytic_Odfmg_Opacity.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 13 17:24:12 2001
 * \brief  Analytic_Odfmg_Opacity test.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "cdi_analytic_test.hh"
#include "c4/ParallelUnitTest.hh"
#include "cdi/CDI.hh"
#include "cdi_analytic/Pseudo_Line_Analytic_Odfmg_Opacity.hh"
#include "cdi_analytic/nGray_Analytic_Odfmg_Opacity.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "parser/Constant_Expression.hh"
#include <sstream>

using namespace std;
using namespace rtt_cdi_analytic;
using namespace rtt_dsxx;
using rtt_cdi::CDI;
using rtt_cdi::OdfmgOpacity;
using rtt_parser::Constant_Expression;
using rtt_parser::Expression;

//---------------------------------------------------------------------------//

bool checkOpacityEquivalence(vector<vector<double>> sigma, vector<double> ref) {
  bool itPasses = true;

  if (sigma.size() != ref.size()) {
    cout << "Mismatch in number of groups: reference " << ref.size()
         << ", from opacity " << sigma.size() << endl;

    return false;
  }

  for (size_t group = 0; group < sigma.size(); group++) {
    for (size_t band = 0; band < sigma[group].size(); band++) {
      if (!soft_equiv(sigma[group][band], ref[group])) {
        itPasses = false;
        cout << "Mismatch in opacities for group " << group << " band " << band
             << endl;
      }
    }
  }
  return itPasses;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void odfmg_test(UnitTest &ut) {
  // group structure
  vector<double> groups = {0.05, 0.5, 5.0, 50.0};

  // band strucutre
  vector<double> bands = {0.0, 0.75, 1.0};

  vector<std::shared_ptr<Analytic_Opacity_Model>> models(3);

  // make a Marshak (user-defined) model for the first group
  models[0].reset(new rtt_cdi_analytic_test::Marshak_Model(100.0));

  // make a Polynomial model for the second group
  models[1].reset(new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
      1.5, 0.0, 0.0, 0.0));

  // make a Constant model for the third group
  models[2].reset(new rtt_cdi_analytic::Constant_Analytic_Opacity_Model(3.0));

  // make an analytic multigroup opacity object for absorption
  nGray_Analytic_Odfmg_Opacity opacity(groups, bands, models,
                                       rtt_cdi::ABSORPTION);

  // check the interface to multigroup opacity
  {
    string desc = "Analytic Odfmg Absorption";

    if (opacity.getOpacityModelType() != rtt_cdi::ANALYTIC_TYPE)
      ITFAILS;
    if (opacity.data_in_tabular_form())
      ITFAILS;
    if (opacity.getReactionType() != rtt_cdi::ABSORPTION)
      ITFAILS;
    if (opacity.getModelType() != rtt_cdi::ANALYTIC)
      ITFAILS;
    if (opacity.getNumTemperatures() != 0)
      ITFAILS;
    if (opacity.getNumDensities() != 0)
      ITFAILS;
    if (opacity.getTemperatureGrid() != vector<double>())
      ITFAILS;
    if (opacity.getDensityGrid() != vector<double>())
      ITFAILS;
    if (opacity.getNumGroups() != 3)
      ITFAILS;
    if (opacity.getNumGroupBoundaries() != 4)
      ITFAILS;
    if (opacity.getNumBands() != 2)
      ITFAILS;
    if (opacity.getNumBandBoundaries() != 3)
      ITFAILS;
    if (opacity.getEnergyPolicyDescriptor() != "odfmg")
      ITFAILS;
    if (opacity.getDataDescriptor() != desc)
      ITFAILS;
    if (opacity.getDataFilename() != string())
      ITFAILS;

    if (opacity.get_Analytic_Model(1, 1) != models[0])
      ITFAILS;
    if (opacity.get_Analytic_Model(2, 1) != models[1])
      ITFAILS;
    if (opacity.get_Analytic_Model(3, 1) != models[2])
      ITFAILS;
  }

  {
    // make an analytic multigroup opacity object for scattering
    nGray_Analytic_Odfmg_Opacity opac(groups, bands, models,
                                      rtt_cdi::SCATTERING);
    string desc = "Analytic Odfmg Scattering";

    if (opac.getDataDescriptor() != desc)
      ITFAILS;
  }
  {
    // make an analytic multigroup opacity object for scattering
    nGray_Analytic_Odfmg_Opacity opac(groups, bands, models, rtt_cdi::TOTAL);
    string desc = "Analytic Odfmg Total";

    if (opac.getDataDescriptor() != desc)
      ITFAILS;
  }

  // check the group structure
  vector<double> get_groups = opacity.getGroupBoundaries();

  if (soft_equiv(get_groups.begin(), get_groups.end(), groups.begin(),
                 groups.end())) {
    ut.passes("Group boundaries match.");
  } else {
    ut.failure("Group boundaries do not match.");
  }

  // check the band structure
  vector<double> get_bands = opacity.getBandBoundaries();

  if (soft_equiv(get_bands.begin(), get_bands.end(), bands.begin(),
                 bands.end())) {
    ut.passes("Band boundaries match.");
  } else {
    ut.failure("Band boundaries do not match.");
  }

  // >>> get opacities

  // scalar density and temperature
  vector<double> ref = {100.0 / 8.0, 1.5, 3.0};

  // load groups * bands opacities; all bands inside each group should be
  // the same
  vector<vector<double>> sigma = opacity.getOpacity(2.0, 3.0);

  // check for each band and group
  bool itPasses = checkOpacityEquivalence(sigma, ref);

  if (itPasses) {
    ostringstream message;
    message << "Analytic multigroup opacities are correct for "
            << "scalar temperature and scalar density.";
    ut.passes(message.str());
  } else {
    ostringstream message;
    message << "Analytic multigroup opacities are NOT correct for "
            << "scalar temperature and scalar density.";
    ut.failure(message.str());
  }

  // scalar density/temperature + vector density/temperature
  vector<double> data_field(3, 2.0);
  vector<vector<vector<double>>> sig_t = opacity.getOpacity(data_field, 3.0);
  vector<vector<vector<double>>> sig_rho = opacity.getOpacity(2.0, data_field);

  itPasses = true;
  for (int i = 0; i < 3; i++) {
    itPasses = itPasses && checkOpacityEquivalence(sig_t[i], ref);
  }

  if (itPasses) {
    ostringstream message;
    message << "Analytic multigroup opacities are correct for "
            << "temperature field and scalar density.";
    ut.passes(message.str());
  } else {
    ostringstream message;
    message << "Analytic multigroup opacities are NOT correct for "
            << "temperature field and scalar density.";
    ut.failure(message.str());
  }

  itPasses = true;
  for (int i = 0; i < 3; i++) {
    itPasses = itPasses && checkOpacityEquivalence(sig_rho[i], ref);
  }

  if (itPasses) {
    ostringstream message;
    message << "Analytic multigroup opacities are correct for "
            << "density field and scalar temperature.";
    ut.passes(message.str());
  } else {
    ostringstream message;
    message << "Analytic multigroup opacities are NOT correct for "
            << "density field and scalar temperature.";
    ut.failure(message.str());
  }

  // Test the get_Analytic_Model() member function.
  std::shared_ptr<Analytic_Opacity_Model const> my_mg_opacity_model =
      opacity.get_Analytic_Model(1);
  std::shared_ptr<Analytic_Opacity_Model const> expected_model(models[0]);

  if (expected_model == my_mg_opacity_model)
    ut.passes(std::string("get_Analytic_Model() returned the expected") +
              std::string(" MG Opacity model."));
  else
    ut.failure(std::string("get_Analytic_Model() did not return the") +
               std::string(" expected MG Opacity model."));

  return;
}

//---------------------------------------------------------------------------//
void test_CDI(UnitTest &ut) {
  // group structure
  vector<double> groups = {0.05, 0.5, 5.0, 50.0};

  // band strucutre
  vector<double> bands = {0.0, 0.75, 1.0};

  vector<std::shared_ptr<Analytic_Opacity_Model>> models(3);

  // make a Marshak (user-defined) model for the first group
  models[0].reset(new rtt_cdi_analytic_test::Marshak_Model(100.0));

  // make a Polynomial model for the second group
  models[1].reset(new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
      1.5, 0.0, 0.0, 0.0));

  // make a Constant model for the third group
  models[2].reset(new rtt_cdi_analytic::Constant_Analytic_Opacity_Model(3.0));

  // make an analytic multigroup opacity object for absorption
  std::shared_ptr<const OdfmgOpacity> odfmg(new nGray_Analytic_Odfmg_Opacity(
      groups, bands, models, rtt_cdi::ABSORPTION));

  // make a CDI object
  CDI cdi;

  // set the multigroup opacity
  cdi.setOdfmgOpacity(odfmg);

  // check the energy groups from CDI
  vector<double> odfmg_groups = CDI::getFrequencyGroupBoundaries();

  if (soft_equiv(odfmg_groups.begin(), odfmg_groups.end(), groups.begin(),
                 groups.end())) {
    ut.passes("CDI Group boundaries match.");
  } else {
    ut.failure("CDI Group boundaries do not match.");
  }

  // check the energy groups from CDI
  vector<double> odfmg_bands = CDI::getOpacityCdfBandBoundaries();

  if (soft_equiv(odfmg_bands.begin(), odfmg_bands.end(), bands.begin(),
                 bands.end())) {
    ut.passes("CDI band boundaries match.");
  } else {
    ut.failure("CDI band boundaries do not match.");
  }

  // do a quick access test for getOpacity

  // scalar density and temperature
  vector<vector<double>> sigma =
      cdi.odfmg(rtt_cdi::ANALYTIC, rtt_cdi::ABSORPTION)->getOpacity(2.0, 3.0);

  vector<double> ref = {100.0 / 8.0, 1.5, 3.0};

  if (checkOpacityEquivalence(sigma, ref)) {
    ostringstream message;
    message << "CDI odfmg opacities are correct for "
            << "scalar temperature and scalar density.";
    ut.passes(message.str());
  } else {
    ostringstream message;
    message << "CDI odfmg opacities are NOT correct for "
            << "scalar temperature and scalar density.";
    ut.failure(message.str());
  }
}

//---------------------------------------------------------------------------//
void packing_test(UnitTest &ut) {
  vector<char> packed;

  // group structure
  vector<double> groups = {0.05, 0.5, 5.0, 50.0};

  // band strucutre
  vector<double> bands = {0.0, 0.75, 1.0};

  {
    vector<std::shared_ptr<Analytic_Opacity_Model>> models(3);

    // make a Polynomial model for the first group
    models[0].reset(new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
        0.0, 100.0, -3.0, 0.0));

    // make a Polynomial model for the second group
    models[1].reset(new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
        1.5, 0.0, 0.0, 0.0));

    // make a Constant model for the third group
    models[2].reset(new rtt_cdi_analytic::Constant_Analytic_Opacity_Model(3.0));

    // make an analytic multigroup opacity object for absorption
    std::shared_ptr<const OdfmgOpacity> odfmg(new nGray_Analytic_Odfmg_Opacity(
        groups, bands, models, rtt_cdi::ABSORPTION));

    // pack it
    packed = odfmg->pack();
  }

  // now unpack it
  nGray_Analytic_Odfmg_Opacity opacity(packed);

  // now check it

  // check the interface to multigroup opacity
  {
    string desc = "Analytic Odfmg Absorption";

    if (opacity.data_in_tabular_form())
      ITFAILS;
    if (opacity.getReactionType() != rtt_cdi::ABSORPTION)
      ITFAILS;
    if (opacity.getModelType() != rtt_cdi::ANALYTIC)
      ITFAILS;
    if (opacity.getNumTemperatures() != 0)
      ITFAILS;
    if (opacity.getNumDensities() != 0)
      ITFAILS;
    if (opacity.getTemperatureGrid() != vector<double>())
      ITFAILS;
    if (opacity.getDensityGrid() != vector<double>())
      ITFAILS;
    if (opacity.getNumGroups() != 3)
      ITFAILS;
    if (opacity.getNumGroupBoundaries() != 4)
      ITFAILS;
    if (opacity.getNumBands() != 2)
      ITFAILS;
    if (opacity.getNumBandBoundaries() != 3)
      ITFAILS;
    if (opacity.getEnergyPolicyDescriptor() != "odfmg")
      ITFAILS;
    if (opacity.getDataDescriptor() != desc)
      ITFAILS;
    if (opacity.getDataFilename() != string())
      ITFAILS;
  }

  // check the group structure
  vector<double> mg_groups = opacity.getGroupBoundaries();

  if (soft_equiv(mg_groups.begin(), mg_groups.end(), groups.begin(),
                 groups.end())) {
    ut.passes("Group boundaries for unpacked ODFMG opacity match.");
  } else {
    ut.failure("Group boundaries for unpacked ODFMG do not match.");
  }

  // check the band structure
  vector<double> get_bands = opacity.getBandBoundaries();

  if (soft_equiv(get_bands.begin(), get_bands.end(), bands.begin(),
                 bands.end())) {
    ut.passes("Band boundaries match.");
  } else {
    ut.failure("Band boundaries do not match.");
  }

  // >>> get opacities

  // scalar density and temperature
  vector<vector<double>> sigma = opacity.getOpacity(2.0, 3.0);
  vector<double> ref = {100.0 / 8.0, 1.5, 3.0};

  if (checkOpacityEquivalence(sigma, ref)) {
    ostringstream message;
    message << "Analytic multigroup opacities for unpacked MG opacity "
            << "are correct for "
            << "scalar temperature and scalar density.";
    ut.passes(message.str());
  } else {
    ostringstream message;
    message << "Analytic multigroup opacities for unpacked MG opacity "
            << "are NOT correct for "
            << "scalar temperature and scalar density.";
    ut.failure(message.str());
  }

  // make sure we catch an assertion showing that we cannot unpack an
  // unregistered opacity
  {
    vector<std::shared_ptr<Analytic_Opacity_Model>> models(3);

    // make a Marshak (user-defined) model for the first group
    models[0].reset(new rtt_cdi_analytic_test::Marshak_Model(100.0));

    // make a Polynomial model for the second group
    models[1].reset(new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
        1.5, 0.0, 0.0, 0.0));

    // make a Constant model for the third group
    models[2].reset(new rtt_cdi_analytic::Constant_Analytic_Opacity_Model(3.0));

    // make an analytic multigroup opacity object for absorption
    std::shared_ptr<const OdfmgOpacity> odfmg(new nGray_Analytic_Odfmg_Opacity(
        groups, bands, models, rtt_cdi::ABSORPTION));

    packed = odfmg->pack();
  }

  // we should catch an assertion when unpacking this because the
  // Marshak_Model is not registered in rtt_cdi::Opacity_Models
  bool caught = false;
  try {
    nGray_Analytic_Odfmg_Opacity nmg(packed);
  } catch (const rtt_dsxx::assertion &err) {
    caught = true;
    ostringstream message;
    message << "Caught the following assertion, " << err.what();
    ut.passes(message.str());
  }
  if (!caught) {
    ut.failure("Failed to catch unregistered analyic model assertion");
  }
}

//---------------------------------------------------------------------------//
void pseudo_line_opacity_test(UnitTest &ut) {
  // group structure
  vector<double> groups = {0.05, 0.5, 5.0, 10.0};

  // band strucutre
  vector<double> bands = {0.0, 0.75, 1.0};

  size_t const number_of_energy_groups = groups.size() - 1;
  size_t const bands_per_group = bands.size() - 1;

  // continuum
  std::shared_ptr<Expression const> const continuum(
      new Constant_Expression(1, 1.0));

  Pseudo_Line_Analytic_Odfmg_Opacity model(groups, bands, rtt_cdi::ABSORPTION,
                                           continuum,
                                           100,   // lines
                                           100.0, // line peak
                                           0.001, // line width
                                           10,    // edges
                                           10,    // edge ratio
                                           1.0,   // Tref
                                           0.0,   // Tpow: no temp dependence
                                           0.0,   // emin keV
                                           10.0,  // emax keV
                                           Pseudo_Line_Base::ROSSELAND,
                                           2000, // quadrature points
                                           1);   // random seed

  vector<vector<double>> opacity = model.getOpacity(1.0, 1.0);

  if (opacity.size() + 1 == groups.size())
    ut.passes("Correct number of groups in pseudo line model");
  else
    ut.failure("NOT correct number of groups in pseudo line model");

  for (unsigned i = 0; i < number_of_energy_groups; ++i) {
    if (opacity[i].size() != bands_per_group) {
      ut.failure("NOT correct number of bands per group "
                 "in pseudo line model");
    }
  }

  vector<vector<vector<double>>> opacities =
      model.getOpacity(vector<double>(2, 1.0), 1.0);

  if (opacities.size() != 2) {
    ut.failure("NOT correct number of opacities in vector T opacity call");
  }

  opacities = model.getOpacity(1.0, vector<double>(2, 1.0));

  if (opacities.size() != 2) {
    ut.failure(std::string("NOT correct number of opacities in vector") +
               std::string(" rho opacity call"));
  }

  if (model.getDataDescriptor() != "Pseudo Line Odfmg Absorption") {
    ut.failure("NOT correct data descriptor");
  }

  // Try pack

  // kgbudge: Doesn't work yet, because we haven't implemented packing for
  // expression trees yet.

  // vector<char> data = model.pack();
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    odfmg_test(ut);
    test_CDI(ut);
    packing_test(ut);
    pseudo_line_opacity_test(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstAnalytic_Odfmg_Opacity.cc
//---------------------------------------------------------------------------//
