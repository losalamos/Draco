//-----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton2/test/tCompton2.cc
 * \author Andrew Till
 * \date   14 Oct 2020
 * \brief  Implementation file for tCompton2
 * \note   Copyright (C) 2017-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "compton2/Compton2.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace rtt_compton2_test {

using rtt_dsxx::soft_equiv;

//----------------------------------------------------------------------------//
// TESTS
//----------------------------------------------------------------------------//

//!  Simple test
void test(rtt_dsxx::UnitTest &ut) {
  // Make true if golds need updating
  const bool do_print = true;

  // Tolerance used for checks
  // TODO: Figure out why I can't use 1e-11 for the tol (outscat)
  const double tol = 2e-6;

  // Start the test.

  std::cout << "\n---------------------------------------------------------\n"
            << "             Test Draco direct CSK routines\n"
            << "---------------------------------------------------------\n";

  // open a small mg opacity file:
  const std::string filename = ut.getTestSourcePath() + "dummy_data_b";
  std::cout << "Attempting to construct a Compton2 object...\n" << std::endl;
  std::unique_ptr<rtt_compton2::Compton2> compton_test;

  try {
    compton_test.reset(new rtt_compton2::Compton2(filename));
  } catch (int /*asrt*/) {
    FAILMSG("Failed to construct a Compton2 object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "\n(...Success!)" << std::endl;

  // Check some of the data in the CSK_generator-opened file:
  const std::vector<double> grp_bds = compton_test->get_Egs();
  const std::vector<double> T_evals = compton_test->get_Ts();

  // Unitless (divided by mec2)
  // NB: These values can be read directly from the ASCII data files
  // (3rd line)
  std::vector<double> grp_bds_gold = {1.57311251e-06, 3.14622503e-04, 7.86556258e-04,
                                      1.57311251e-03, 3.14622503e-02};
  // (scattered throughout data file)
  std::vector<double> T_evals_gold = {1.57311251e-05, 1.57311251e-04, 3.30353629e-04,
                                      6.60707256e-04};

  // Sizes
  const size_t num_groups_gold = grp_bds_gold.size() - 1U;
  const size_t num_T_evals_gold = T_evals_gold.size();
  const size_t num_evals_gold = 3; // in_lin, out_lin, diff_nonlin
  const size_t num_leg_moments_gold = 2;
  // A point is a (Legendre moment, evaluation) pair
  // first eval (in_lin) has all Leg moments and all others have only the 0th moment
  const size_t num_points_gold = num_leg_moments_gold + (num_evals_gold - 1U);

  // Multiply by electron rest-mass energy (keV; using CSK value)
  const double mec2 = 510.998;
  for (size_t i = 0; i < grp_bds_gold.size(); ++i) {
    grp_bds_gold[i] *= mec2;
  }
  for (size_t i = 0; i < T_evals_gold.size(); ++i) {
    T_evals_gold[i] *= mec2;
  }

  ut.check(grp_bds.size() == (num_groups_gold + 1U), "checked size of group bounds vector");
  ut.check(std::equal(grp_bds.begin(), grp_bds.end(), grp_bds_gold.begin(),
                      [tol](double a, double b) -> bool { return soft_equiv(a, b, tol); }),
           "checked group boundaries");

  ut.check(T_evals.size() == num_T_evals_gold, "checked size of temperature evals vector");
  ut.check(std::equal(T_evals.begin(), T_evals.end(), T_evals_gold.begin(),
                      [tol](double a, double b) -> bool { return soft_equiv(a, b, tol); }),
           "checked temperature grid");

  // Test size accessor functions
  ut.check(compton_test->get_num_temperatures() == num_T_evals_gold,
           "checked number of temperatures");
  ut.check(compton_test->get_num_groups() == num_groups_gold, "checked number of groups");
  ut.check(compton_test->get_num_leg_moments() == num_leg_moments_gold,
           "checked number of Legendre moments");
  ut.check(compton_test->get_num_evals() == num_evals_gold, "checked number of evaluations");
  ut.check(compton_test->get_num_points() == num_points_gold, "checked number of points");
  ut.check(compton_test->get_highest_leg_moment() == (num_leg_moments_gold - 1U),
           "checked highest Legendre moment");

  if (ut.numFails == 0) {
    std::cout << "\nCorrectly read sizes, group bounds, and electron temps!" << std::endl;
  }

  // Test data retrieval: interpolate to a grid point in temperature
  {
    const double interp_T_keV = T_evals[num_T_evals_gold - 1U];
    std::cout << "Testing interpolation at T = " << interp_T_keV << " keV\n";

    const size_t G = compton_test->get_num_groups();
    const size_t L = compton_test->get_num_leg_moments();
    std::vector<double> inscat(G * G * L, -1.0);
    std::vector<double> outscat(G, -1.0);
    std::vector<double> nl_diff(G, -1.0);

    // Returns flattened 3D inscat array with order [moment, group_to, group_from]
    compton_test->interp_dense_inscat(inscat, interp_T_keV, G);

    // Returns 1D outscat array [group_from]
    compton_test->interp_linear_outscat(outscat, interp_T_keV);

    // Use dummy flux for nonlinear component
    const double phival = 2.0;
    const double phiscale = 0.1;
    std::vector<double> phi(compton_test->get_num_groups(), phival);
    // Returns 1D nonlinear difference array [group_from] (performs matrix multiplication with phi)
    compton_test->interp_nonlin_diff_and_add(nl_diff, interp_T_keV, phi, phiscale);

    // Print result (useful if golds need updating)
    if (do_print) {
      [](const std::vector<double> &v) {
        auto print = [](double a) { std::cout << a << ", "; };
        std::cout << std::setprecision(14);
        std::cout << "\ninscat\n";
        std::for_each(v.begin(), v.end(), print);
        std::cout << std::endl;
      }(inscat);

      [](const std::vector<double> &v) {
        auto print = [](double a) { std::cout << a << ", "; };
        std::cout << std::setprecision(14);
        std::cout << "\noutscat\n";
        std::for_each(v.begin(), v.end(), print);
        std::cout << std::endl;
      }(outscat);

      [](const std::vector<double> &v) {
        auto print = [](double a) { std::cout << a << ", "; };
        std::cout << std::setprecision(14);
        std::cout << "\nnl_diff\n";
        std::for_each(v.begin(), v.end(), print);
        std::cout << std::endl;
      }(nl_diff);
    }

    // NB: These dense values come directly from the ASCII data files
    // flattened 1D array with ordering (slowest) [eval, gfrom, gto, leg] (fastest)
    std::vector<double> raw_gold = {
        // out_lin
        1.023686968316, 0.003573655955675, 0.1408569058113, -0.003385368827333, 0, 0, 0, 0,
        0.01368733137026, -0.000337289996441, 0.8580132416506, 0.00346237220459, 0.1150544176069,
        -0.002695776219916, 0, 0, 0, 0, 0.03016927525143, -0.000745487878607, 0.7395544180756,
        0.003537535339196, 0.08790804121747, -0.001974567477495, 0, 0, 0, 0, 0.02760683687763,
        -0.0006747774819197, 0.7154380590368, 0.002500085526454,
        // in_lin
        1.029549163582, 0.003360819879047, 0.1682617986774, -0.004498320348364, 0, 0, 0, 0,
        0.01155720790042, -0.000254680100335, 0.8615393494377, 0.003337305544177, 0.1370933267956,
        -0.003577889661548, 0, 0, 0, 0, 0.02546236651279, -0.0005628307049739, 0.74122185377,
        0.003480918190072, 0.104392775669, -0.002618892286625, 0, 0, 0, 0, 0.02330836203736,
        -0.0005092109709329, 0.7158405069922, 0.002485525111488,
        // out_nonlin
        2.417350052736e-20, 5.333857003773e-23, 1.569215319107e-22, -3.474039416897e-24, 0, 0, 0, 0,
        2.025140617993e-22, -5.27716088378e-24, 4.991731616824e-22, 1.576887628939e-24,
        7.937214293507e-24, -1.661082680445e-25, 0, 0, 0, 0, 1.15611911544e-23, -3.08379217909e-25,
        2.906510852958e-23, 1.161410570085e-25, 3.652307822369e-25, -6.9213492314e-27, 0, 0, 0, 0,
        6.471151078044e-25, -1.788886184269e-26, 9.067645193159e-25, 3.288019868804e-27,
        // in_nonlin
        2.412582001649e-20, 5.506360234252e-23, 1.858984933376e-22, -4.598103855448e-24, 0, 0, 0, 0,
        1.70028840385e-22, -3.978450052949e-24, 4.964897275093e-22, 1.671777660171e-24,
        9.356135881693e-24, -2.194753402066e-25, 0, 0, 0, 0, 9.680369923644e-24, -2.32267253788e-25,
        2.875946793361e-23, 1.265594484598e-25, 4.274264612141e-25, -9.142968009692e-27, 0, 0, 0, 0,
        5.392751048597e-25, -1.344066881751e-26, 8.87128629256e-25, 3.986007245261e-27};

    // Slice up raw gold and rescale
    const double mtocm = 100.0;
    const double classical_electron_radius = mtocm * rtt_units::classicalElectronRadiusSI; // cm
    const double valscale = 0.25 * 2 * rtt_units::PI * classical_electron_radius *
                            classical_electron_radius * rtt_units::AVOGADRO;
    std::cout << "Scaling raw opacities by " << valscale << '\n';

    // Fill inscat_gold (reorder)
    // TODO: FIX REORDERING
    std::vector<double> inscat_gold(G * G * L, 0.0);
    const size_t in_offset = G * G * L;
    for (size_t gfrom = 0; gfrom < G; ++gfrom) {
      for (size_t gto = 0; gto < G; ++gto) {
        for (size_t m = 0; m < L; ++m) {
          const size_t i_raw = in_offset + m + L * (gto + G * gfrom);
          const size_t i = gfrom + G * (gto + G * m);
          inscat_gold[i] = valscale * raw_gold[i_raw];
        }
      }
    }

    // Fill outscat_gold (reorder and sum)
    std::vector<double> outscat_gold(G, 0.0);
    const size_t out_offset = 0;
    const size_t P0 = 0;
    for (size_t gfrom = 0; gfrom < G; ++gfrom) {
      double sum = 0.0;
      for (size_t gto = 0; gto < G; ++gto) {
        const size_t i_raw = out_offset + P0 + L * (gto + G * gfrom);
        sum += raw_gold[i_raw];
      }
      outscat_gold[gfrom] = valscale * sum;
    }

    std::vector<double> nl_diff_gold(G, 0.0);

    // Print result (useful if golds need updating)
    if (do_print) {
      [](const std::vector<double> &v) {
        auto print = [](double a) { std::cout << a << ", "; };
        std::cout << std::setprecision(14);
        std::cout << "\ninscat_gold\n";
        std::for_each(v.begin(), v.end(), print);
        std::cout << std::endl;
      }(inscat_gold);

      [](const std::vector<double> &v) {
        auto print = [](double a) { std::cout << a << ", "; };
        std::cout << std::setprecision(14);
        std::cout << "\noutscat_gold\n";
        std::for_each(v.begin(), v.end(), print);
        std::cout << std::endl;
      }(outscat_gold);

      [](const std::vector<double> &v) {
        auto print = [](double a) { std::cout << a << ", "; };
        std::cout << std::setprecision(14);
        std::cout << "\nnl_diff_gold\n";
        std::for_each(v.begin(), v.end(), print);
        std::cout << std::endl;
      }(nl_diff_gold);
    }

    // Print diff (useful if golds need updating)
    if (do_print) {
      [tol](const std::vector<double> &v, const std::vector<double> &r) {
        std::cout << std::setprecision(14);
        std::cout << "\ninscat / inscat_gold - 1\n";
        for (size_t i = 0; i < v.size(); ++i)
          std::cout << (v[i] - r[i]) / (std::fabs(r[i]) + tol) << ", ";
        std::cout << std::endl;
      }(inscat, inscat_gold);

      [tol](const std::vector<double> &v, const std::vector<double> &r) {
        std::cout << std::setprecision(14);
        std::cout << "\noutscat / outscat_gold - 1\n";
        for (size_t i = 0; i < v.size(); ++i)
          std::cout << (v[i] - r[i]) / (std::fabs(r[i]) + tol) << ", ";
        std::cout << std::endl;
      }(outscat, outscat_gold);
    }
    // TODO: Hack print of values from cskrw.cc
    // TODO: Check sparsity structure of interp or data
    // TODO: Check temperature interpolation coeff (verified correct)

    // TODO: Figure out why inscat has unexpected non-zeros
    ut.check(std::equal(inscat.begin(), inscat.end(), inscat_gold.begin(),
                        [tol](double a, double b) -> bool { return soft_equiv(a, b, tol); }),
             "checked data retrieval for inscat");

    // TODO: Figure out why outscat is slightly off
    ut.check(std::equal(outscat.begin(), outscat.end(), outscat_gold.begin(),
                        [tol](double a, double b) -> bool { return soft_equiv(a, b, tol); }),
             "checked data retrieval for outscat");

    // Enough operations are done on the nl_diff in the converter
    // that using the raw values is not feasible
    // TODO: Hack print of nl_diff and use for gold
  }
}

} // namespace rtt_compton2_test

//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    rtt_compton2_test::test(ut);
  }
  UT_EPILOG(ut);
}

//----------------------------------------------------------------------------//
// End of test/tCompton2.cc
//----------------------------------------------------------------------------//
