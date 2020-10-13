//-----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/test/t2Compton.cc
 * \author Andrew Till
 * \date   2020 Oct 14
 * \brief  Implementation file for t2Compton
 * \note   Copyright (C) 2017-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "compton/Compton.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace rtt_compton_test {

using rtt_dsxx::soft_equiv;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//

//!  Tests the Compton constructor and a couple of access routines.
void compton_file_test(rtt_dsxx::UnitTest &ut) {

  const bool do_print = true;

  // Tolerance used for checks
  const double tol = 1e-11;

  // Start the test.

  std::cout << "\n---------------------------------------------------------\n"
            << "   Test Draco code calling CSK_generator routines\n"
            << "---------------------------------------------------------\n";

  // open a small mg opacity file:
  const std::string filename = ut.getTestSourcePath() + "../../compton2/test/dummy_data";
  std::cout << "Attempting to construct a Compton object...\n" << std::endl;
  std::unique_ptr<rtt_compton::Compton> compton_test;

  try {
    compton_test.reset(new rtt_compton::Compton(filename));
  } catch (int /*asrt*/) {
    FAILMSG("Failed to construct a Compton object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "\n(...Success!)" << std::endl;

  // Check some of the data in the CSK_generator-opened file:
  const std::vector<double> grp_bds = compton_test->get_group_bounds();
  const std::vector<double> T_evals = compton_test->get_etemp_pts();

  // Unitless
  std::vector<double> grp_bds_gold = {1.57311251e-06, 3.14622503e-04, 7.86556258e-04,
                                      1.57311251e-03, 3.14622503e-02};
  std::vector<double> T_evals_gold = {1.57311251e-05, 1.57311251e-04, 3.30353629e-04,
                                      6.60707256e-04};
  // First and last temperature from line 2 of the csk input file
  std::vector<double> line2_Ts_gold = {1.41580126e-05, 7.26777982e-04};

  // Multiply by electron rest-mass energy (keV; using CSK value)
  const double mec2 = 510.998;
  for (size_t i = 0; i < grp_bds_gold.size(); ++i) {
    grp_bds_gold[i] *= mec2;
  }
  // Interface does not scale temperatures, so no multiplication needed:
  /*
  for (size_t i = 0; i < T_evals_gold.size(); ++i)
  {
      T_evals_gold[i] *= mec2;
  }
  line2_Ts_gold[0] *= mec2;
  line2_Ts_gold[1] *= mec2;
  */

  Ensure(grp_bds.size() == grp_bds_gold.size());
  ut.check(std::equal(grp_bds.begin(), grp_bds.end(), grp_bds_gold.begin(),
                      [tol](double a, double b) -> bool { return soft_equiv(a, b, tol); }),
           "checked group boundaries");

  Ensure(T_evals.size() == T_evals_gold.size());
  ut.check(std::equal(T_evals.begin(), T_evals.end(), T_evals_gold.begin(),
                      [tol](double a, double b) -> bool { return soft_equiv(a, b, tol); }),
           "checked temperature grid");

  if (!soft_equiv(compton_test->get_min_etemp(), line2_Ts_gold[0]))
    FAILMSG("Min etemp read incorrectly!");
  if (!soft_equiv(compton_test->get_max_etemp(), line2_Ts_gold[1]))
    FAILMSG("Max etemp read incorrectly!");

  if (ut.numFails == 0) {
    std::cout << "\nCorrectly read group bounds and electron temps!" << std::endl;
  }

  // Interpolate
  const double alpha0 = 0.4;
  const double interp_T0 = alpha0 * T_evals[0] + (1.0 - alpha0) * T_evals[1];
  std::cout << "Testing interpolation at T = " << (interp_T0 * mec2) << " keV\n";

  {
    std::vector<std::vector<std::vector<std::vector<double>>>> interp_data =
        compton_test->interpolate_csk(interp_T0);

    const size_t sz0 = interp_data.size();
    const size_t sz1 = interp_data[0].size();
    const size_t sz2 = interp_data[0][0].size();
    const size_t sz3 = interp_data[0][0][0].size();
    const size_t tot_size = sz0 * sz1 * sz2 * sz3;
    std::vector<double> flat_interp_data(tot_size, 0.0);

    for (size_t i0 = 0; i0 < sz0; ++i0) {
      for (size_t i1 = 0; i1 < sz1; ++i1) {
        for (size_t i2 = 0; i2 < sz2; ++i2) {
          for (size_t i3 = 0; i3 < sz3; ++i3) {
            const size_t i = i3 + sz3 * (i2 + sz2 * (i1 + sz1 * (i0)));
            flat_interp_data[i] = interp_data[i0][i1][i2][i3];
          }
        }
      }
    }

    // Test sizes
    const size_t num_evals_gold = 4; // in_lin, out_lin, in_nonlin, out_nonlin
    const size_t num_groups_gold = grp_bds_gold.size() - 1U;
    const size_t num_leg_moments_gold = 2;
    ut.check(sz0 == num_evals_gold, "tested evals size");
    ut.check((sz1 == num_groups_gold) && (sz2 == num_groups_gold), "tested groups size");
    ut.check(sz3 == num_leg_moments_gold, "tested Legendre moments size");

    // Print result (useful if golds need updating)
    if (do_print) {
      [](const std::vector<double> &v) {
        auto print = [](double a) { std::cout << a << ", "; };
        std::cout << std::setprecision(14);
        std::cout << '\n';
        std::for_each(v.begin(), v.end(), print);
        std::cout << std::endl;
      }(flat_interp_data);
    }

    // TODO: Test data and do another interp in center

    std::vector<double> flat_interp_gold = {1.4895048336024,
                                            0.0013157633414251,
                                            0.04145112467826,
                                            -0.00099234375153514,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.048935718979052,
                                            -0.001201759773222,
                                            1.1427518955345,
                                            0.0020642528309125,
                                            0.015831053793403,
                                            -0.00036191022639507,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.060076056533251,
                                            -0.0014978191061645,
                                            0.95368017957811,
                                            0.0024195173860433,
                                            0.0042167573598475,
                                            -3.7346606985974e-05,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.12606997604209,
                                            -0.0029756533421552,
                                            0.70308466958259,
                                            0.0057984236887092,
                                            1.4902075219537,
                                            0.0012890343406445,
                                            0.044522269416302,
                                            -0.0011149333854142,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.046790157320168,
                                            -0.0011165953164094,
                                            1.1422078133786,
                                            0.002071967455139,
                                            0.016959526258933,
                                            -0.00040597846461434,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.055452234941061,
                                            -0.0013152449665463,
                                            0.95025760553927,
                                            0.002535497503532,
                                            0.0044495495432538,
                                            -4.5183679054704e-05,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0.11605523138559,
                                            -0.0025883973572953,
                                            0.67868811344021,
                                            0.0063593874505928,
                                            4.6630056595372e-18,
                                            8.5759248352589e-23,
                                            1.1673356563462e-22,
                                            -2.6167370751723e-24,
                                            0,
                                            0,
                                            0,
                                            0,
                                            1.5017646436551e-22,
                                            -3.9459638253177e-24,
                                            9.1490085641636e-21,
                                            5.3409650219456e-23,
                                            5.484626607511e-24,
                                            -1.0596457235234e-25,
                                            0,
                                            0,
                                            0,
                                            0,
                                            1.0314177969686e-23,
                                            -2.9818990024258e-25,
                                            4.2638762788237e-22,
                                            8.5042714161135e-25,
                                            1.8073497542315e-25,
                                            -2.4297847629037e-27,
                                            0,
                                            0,
                                            0,
                                            0,
                                            7.6278018408129e-25,
                                            -2.5980363364632e-26,
                                            1.5794222622867e-23,
                                            1.5887774510949e-25,
                                            4.6615665467427e-18,
                                            1.7873552446889e-22,
                                            1.2663879157434e-22,
                                            -2.9982049815629e-24,
                                            0,
                                            0,
                                            0,
                                            0,
                                            1.3770981199909e-22,
                                            -3.4369299797817e-24,
                                            9.1220263142097e-21,
                                            5.4325270040856e-23,
                                            5.890752806507e-24,
                                            -1.2076068345897e-25,
                                            0,
                                            0,
                                            0,
                                            0,
                                            9.5432783203125e-24,
                                            -2.6517064022462e-25,
                                            4.2104110096026e-22,
                                            1.0546962404339e-24,
                                            1.9308590430364e-25,
                                            -2.8485633209286e-27,
                                            0,
                                            0,
                                            0,
                                            0,
                                            7.0910968932275e-25,
                                            -2.3423156364316e-26,
                                            1.4873126610905e-23,
                                            1.8374675563093e-25};

    ut.check(std::equal(flat_interp_data.begin(), flat_interp_data.end(), flat_interp_gold.begin(),
                        [tol](double a, double b) -> bool { return soft_equiv(a, b, tol); }),
             "checked data interpolation");
  }

  return;

  // try "interpolating" at one of the exact eval points in the test library,
  // and check the result:
  const double test_etemp = 4.87227167e-04;
  std::vector<std::vector<std::vector<double>>> interp_data =
      compton_test->interpolate_csk(test_etemp)[0];
  // get interpolated nu_ratios
  std::vector<std::vector<double>> interp_nu_data = compton_test->interpolate_nu_ratio(test_etemp);

  // Check the size of the returned data:
  Ensure(interp_nu_data.size() == 1);
  Ensure(interp_nu_data[0].size() == 1);
  Ensure(interp_data.size() == 1);
  Ensure(interp_data[0].size() == 1);
  Ensure(interp_data[0][0].size() == 4);

  // Check that the data is actually correct:
  if (!soft_equiv(interp_data[0][0][0], 4.45668383e+00))
    ITFAILS;
  if (!soft_equiv(interp_data[0][0][1], 3.17337784e-01))
    ITFAILS;
  if (!soft_equiv(interp_data[0][0][2], 4.50133379e-01))
    ITFAILS;
  if (!soft_equiv(interp_data[0][0][3], 3.59663442e-02))
    ITFAILS;

  if (!soft_equiv(interp_nu_data[0][0], 1.5000000e+00))
    ITFAILS;

  const double test_etemp2 = 6.75507064e-04;
  // get interpolated nu_ratio for a different
  std::vector<std::vector<double>> interp_nu_data2 =
      compton_test->interpolate_nu_ratio(test_etemp2);

  if (!soft_equiv(interp_nu_data2[0][0], 2.0000000e+00))
    ITFAILS;

  // get the number of xi evals in the library (we know it should be 4)
  if (compton_test->get_num_xi() != 4)
    ITFAILS;

  if (compton_test->get_num_groups() != 1)
    ITFAILS;

  if (ut.numFails == 0)
    std::cout << "\nCorrectly read multigroup data points!" << std::endl;

  if (ut.numFails == 0) {
    PASSMSG("Successfully linked Draco against CSK_generator.");
  } else {
    FAILMSG("Did not successfully link Draco against CSK_generator.");
  }
}

//------------------------------------------------------------------------------------------------//
//!  Tests Compton's error-handling on a non-existent file.
void compton_fail_test(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n---------------------------------------------------------\n"
            << "    Test Compton bad file handling    \n"
            << "---------------------------------------------------------\n";
  // open a small mg opacity file:
  std::string filename = ut.getTestSourcePath() + "non_existent.compton";
  std::cout << "Testing with a non-existent file...\n" << std::endl;
  std::unique_ptr<rtt_compton::Compton> compton_test;

  bool caught = false;
  try {
    compton_test.reset(new rtt_compton::Compton(filename));
  } catch (rtt_dsxx::assertion &asrt) {
    std::cout << "Draco exception thrown: " << asrt.what() << std::endl;
    // We successfully caught the bad file!
    caught = true;
  } catch (const int & /*asrt*/) {
    std::cout << "CSK exception thrown. " << std::endl;
    // We successfully caught the bad file!
    caught = true;
  }

  if (!caught)
    ITFAILS;

  if (ut.numFails == 0) {
    PASSMSG("Successfully caught a CSK_generator exception.");
  } else {
    FAILMSG("Did not successfully catch a CSK_generator exception.");
  }
}

//------------------------------------------------------------------------------------------------//
//!  Tests Comptohn interface to LLNL-style Compton data.
void llnl_compton_test(rtt_dsxx::UnitTest &ut) {
  // Start the test.s
  std::cout << "\n---------------------------------------------------------\n"
            << " Test Draco calling CSK_generator LLNL-style routines \n"
            << "---------------------------------------------------------\n";

  const std::string filename = ut.getTestSourcePath() + "llnl_ascii.compton";
  const bool llnl_style = true;
  std::cout << "Attempting to construct a Compton object...\n" << std::endl;
  std::unique_ptr<rtt_compton::Compton> compton_test;
  try {
    compton_test.reset(new rtt_compton::Compton(filename, llnl_style));
  } catch (int /*asrt*/) {
    FAILMSG("Failed to construct an LLNL-style Compton object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "\n(...Success!)" << std::endl;

  // Test the two sets of interpolators, and make sure they match.

  // Use temperatures that are very close to the library eval points:
  const std::vector<double> cell_temps = {1.50001, 2.49999};
  // Ensure densities are unit; these are only folded into the evaluations when
  // the data is pre-interpolated (not in the on-the-fly case)
  const std::vector<double> cell_dens = {1.0, 1.0};

  compton_test->interpolate_precycle(cell_temps, cell_dens);

  const std::vector<double> test_freq = {12.4233, 183.43};

  // reference sol'n for erec:
  const std::vector<double> ref_sigc = {1.91083e-01, 1.25551e-01, 1.91042e-01, 1.25395e-01};
  const std::vector<double> ref_erec = {-1.22843e-02, -2.00713e-01, -5.00821e-03, -1.974413e-01};

  double otf_sigcval;
  double otf_erecval;
  double pre_sigcval;
  double pre_erecval;

  size_t i = 0;
  for (size_t k = 0; k < cell_temps.size(); k++) {
    for (size_t j = 0; j < test_freq.size(); j++) {
      // sigma_compton value:
      otf_sigcval = compton_test->interpolate_sigc(cell_temps[k], test_freq[j]);
      pre_sigcval = compton_test->interpolate_cell_sigc(static_cast<int64_t>(k), test_freq[j]);

      // expected relative energy change value:
      otf_erecval = compton_test->interpolate_erec(cell_temps[k], test_freq[j]);
      pre_erecval = compton_test->interpolate_cell_erec(static_cast<int64_t>(k), test_freq[j]);

      // compare the values to each other for consistency:
      FAIL_IF_NOT(soft_equiv(otf_sigcval, pre_sigcval));
      FAIL_IF_NOT(soft_equiv(otf_erecval, pre_erecval));

      // compare the values to the expected answer for accuracy:
      // (use a loose tolerance, because our points are close -- but not
      // equal to -- the evaluation points in the library)
      FAIL_IF_NOT(soft_equiv(otf_erecval, ref_erec[i], 1.0e-4));
      FAIL_IF_NOT(soft_equiv(otf_sigcval, ref_sigc[i], 1.0e-4));

      i++;
    }
  }

  if (ut.numFails == 0) {
    PASSMSG("Successfully read an LLNL-style CSK library.");
  } else {
    FAILMSG("Did not successfully read an LLNL-style CSK library.");
  }
}
} // namespace rtt_compton_test

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    rtt_compton_test::compton_file_test(ut);
    //rtt_compton_test::compton_fail_test(ut);
    //rtt_compton_test::llnl_compton_test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// End of test/tCompton.cc
//------------------------------------------------------------------------------------------------//
