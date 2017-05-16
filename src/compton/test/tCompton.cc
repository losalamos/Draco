//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton/test/tCompton.cc
 * \author Kendra Keady
 * \date   2017 Feb 10
 * \brief  Implementation file for tCompton
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#include "compton/Compton.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace rtt_compton_test {

using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

//!  Tests the Compton constructor and a couple of access routines.
void compton_file_test(rtt_dsxx::UnitTest &ut) {
  // Start the test.

  std::cout << "\n---------------------------------------------------------"
            << std::endl;
  std::cout << "   Test Draco code calling CSK_generator routines" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;

  // open a small mg opacity file:
  const std::string filename = "mg_ascii.compton";
  std::cout << "Attempting to construct a Compton object...\n" << std::endl;
  std::shared_ptr<rtt_compton::Compton> compton_test;

  try {
    compton_test.reset(new rtt_compton::Compton(filename));
  } catch (int asrt) {
    FAILMSG("Failed to construct a Compton object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "\n(...Success!)" << std::endl;

  // Check some of the data in the CSK_generator-opened file:
  const std::vector<double> grp_bds = compton_test->get_group_bounds();
  const std::vector<double> etemp_evals = compton_test->get_etemp_pts();

  Ensure(grp_bds.size() == 2);
  Ensure(etemp_evals.size() == 7);

  if (!soft_equiv(grp_bds[0], 3.91389000e-02))
    FAILMSG("Lower group bound read incorrectly!");
  if (!soft_equiv(grp_bds[1], 5.87084000e-02))
    FAILMSG("Upper group bound read incorrectly!");

  if (!soft_equiv(etemp_evals[0], 1.76377944e-05))
    FAILMSG("Etemp eval read incorrectly!");
  if (!soft_equiv(etemp_evals[6], 6.75507064e-04))
    FAILMSG("Etemp eval read incorrectly!");

  if (!soft_equiv(compton_test->get_min_etemp(), 1e-6))
    FAILMSG("Min etemp read incorrectly!");
  if (!soft_equiv(compton_test->get_max_etemp(), 7e-4))
    FAILMSG("Max etemp read incorrectly!");

  if (ut.numFails == 0) {
    std::cout << "\nCorrectly read group bounds and electron temps!"
              << std::endl;
  }

  // try "interpolating" at one of the exact eval points in the test library,
  // and check the result:
  const double test_etemp = 4.87227167e-04;
  std::vector<std::vector<std::vector<double>>> interp_data =
      compton_test->interpolate(test_etemp);

  // Check the size of the returned data:
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

void const_compton_file_test(rtt_dsxx::UnitTest &ut) {
  // Start the test.

  std::cout << "\n---------------------------------------------------------"
            << std::endl;
  std::cout << " Test Draco code calling CSK_generator routines -- const "
            << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;

  // open a small mg opacity file:
  const std::string filename = "mg_ascii.compton";
  std::cout << "Attempting to construct a const Compton object...\n"
            << std::endl;
  std::shared_ptr<const rtt_compton::Compton> compton_test;

  try {
    compton_test.reset(new const rtt_compton::Compton(filename));
  } catch (int asrt) {
    FAILMSG("Failed to construct a Compton object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "\n(...Success!)" << std::endl;

  // Check some of the data in the CSK_generator-opened file:
  const std::vector<double> grp_bds = compton_test->get_group_bounds();
  const std::vector<double> etemp_evals = compton_test->get_etemp_pts();

  Ensure(grp_bds.size() == 2);
  Ensure(etemp_evals.size() == 7);

  if (!soft_equiv(grp_bds[0], 3.91389000e-02))
    FAILMSG("Lower group bound read incorrectly!");
  if (!soft_equiv(grp_bds[1], 5.87084000e-02))
    FAILMSG("Upper group bound read incorrectly!");

  if (!soft_equiv(etemp_evals[0], 1.76377944e-05))
    FAILMSG("Etemp eval read incorrectly!");
  if (!soft_equiv(etemp_evals[6], 6.75507064e-04))
    FAILMSG("Etemp eval read incorrectly!");

  if (!soft_equiv(compton_test->get_min_etemp(), 1e-6))
    FAILMSG("Min etemp read incorrectly!");
  if (!soft_equiv(compton_test->get_max_etemp(), 7e-4))
    FAILMSG("Max etemp read incorrectly!");

  if (ut.numFails == 0) {
    std::cout << "\nCorrectly read group bounds and electron temps!"
              << std::endl;
  }

  // try "interpolating" at one of the exact eval points in the test library,
  // and check the result:
  const double test_etemp = 4.87227167e-04;
  std::vector<std::vector<std::vector<double>>> interp_data =
      compton_test->interpolate(test_etemp);

  // Check the size of the returned data:
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

  // get the number of xi evals in the library (we know it should be 4)
  if (compton_test->get_num_xi() != 4)
    ITFAILS;

  if (compton_test->get_num_groups() != 1)
    ITFAILS;

  if (ut.numFails == 0)
    std::cout << "\nCorrectly read multigroup data points!" << std::endl;

  if (ut.numFails == 0) {
    PASSMSG("Successfully linked Draco against CSK_generator (const version).");
  } else {
    FAILMSG("Did not successfully link Draco against CSK_generator (const "
            "version).");
  }
}

//!  Tests the Compton mg build capability
void compton_build_test(rtt_dsxx::UnitTest &ut) {
  // Start the test.

  std::cout << "\n---------------------------------------------------------"
            << std::endl;
  std::cout << "Test Draco call to CSK_generator mg opacity builder"
            << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;

  // open a small pointwise opacity file:
  const std::string filename = "lagrange_csk_ascii.compton";
  std::cout << "Attempting to construct a Compton object..." << std::endl;

  // make an uninitialized pointer...
  std::shared_ptr<rtt_compton::Compton> compton_test;

  // make a small fake group structure to pass in:
  const std::vector<double> test_groups = {20.0, 30.0, 40.0, 50.0, 60.0};

  // set the number of angular points to retrieve (legendre or otherwise)
  const size_t nxi = 3;

  try {
    // (This call has some output of its own, so we print some newlines
    // around it)
    std::cout << "\n\n";
    compton_test.reset(new rtt_compton::Compton(filename, test_groups, nxi));
    std::cout << "\n\n";
  } catch (rtt_dsxx::assertion &asrt) {
    FAILMSG("Failed to construct a Compton object!");
    // if construction fails, there is no reason to continue testing...
    return;
  }
  std::cout << "(...Success!)" << std::endl;

  // Check some of the data in the CSK_generator-formed data file:
  const std::vector<double> grp_bds = compton_test->get_group_bounds();

  if (!soft_equiv(grp_bds.begin(), grp_bds.end(), test_groups.begin(),
                  test_groups.end()))
    ITFAILS;

  if (!soft_equiv(compton_test->get_min_etemp(), 0.0))
    FAILMSG("Min etemp read incorrectly!");
  if (!soft_equiv(compton_test->get_max_etemp(), 1.0))
    FAILMSG("Max etemp read incorrectly!");

  if (ut.numFails == 0) {
    std::cout << "\nCorrectly stored group bounds and electron temps! "
              << std::endl;
  }

  if (ut.numFails == 0) {
    PASSMSG("Successfully built a CSK_generator mg library.");
  } else {
    FAILMSG("Did not successfully build a CSK_generator mg library.");
  }
}

//!  Tests Compton's error-handling on a non-existent file.
void compton_fail_test(rtt_dsxx::UnitTest &ut) {
  std::cout << "\n---------------------------------------------------------"
            << std::endl;
  std::cout << "    Test Compton bad file handling    " << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  // open a small mg opacity file:
  std::string filename = "non_existent.compton";
  std::cout << "Testing with a non-existent file...\n" << std::endl;
  std::shared_ptr<rtt_compton::Compton> compton_test;

  bool caught = false;
  try {
    compton_test.reset(new rtt_compton::Compton(filename));
  } catch (rtt_dsxx::assertion &asrt) {
    std::cout << "Draco exception thrown: " << asrt.what() << std::endl;
    // We successfully caught the bad file!
    caught = true;
  } catch (const int &asrt) {
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
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    rtt_compton_test::compton_file_test(ut);
    rtt_compton_test::const_compton_file_test(ut);
    rtt_compton_test::compton_build_test(ut);
    rtt_compton_test::compton_fail_test(ut);
  }
  UT_EPILOG(ut);
}
