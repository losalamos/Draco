//----------------------------------*-C++-*----------------------------------//
/*! \file   timestep/test/tstTimeStep.cc
 *  \brief  A driver for the time-step manager test facility.
 *  \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *          All rights reserved.  */
//---------------------------------------------------------------------------//

#include "dummy_package.hh"
#include "c4/ParallelUnitTest.hh"
#include "c4/global.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "timestep/field_ts_advisor.hh"
#include "timestep/fixed_ts_advisor.hh"
#include "timestep/ratio_ts_advisor.hh"
#include "timestep/target_ts_advisor.hh"
#include "timestep/ts_manager.hh"
#include <sstream>

// forward declaration
void run_tests(rtt_dsxx::UnitTest &ut);
void check_field_ts_advisor(rtt_dsxx::UnitTest &ut);

//---------------------------------------------------------------------------//
// Main program
//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try { // Run the tests...
    run_tests(ut);
    if (rtt_c4::node() == 0)
      check_field_ts_advisor(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// Actual tests go here.
//---------------------------------------------------------------------------//

void run_tests(rtt_dsxx::UnitTest &ut) {
  using namespace rtt_timestep;
  using rtt_timestep_test::dummy_package;

  // Initial values;
  double graphics_time(10.0);
  double dt_min(0.000001);
  double dt_max(100000.0);
  double override_dt(1.0);
  double dt(1.0);
  double time(0.);

  int icycle_first(1);
  int icycle_last(3);

  bool override_flag(false);

  ts_manager mngr;
  dummy_package xxx(mngr);

  // Set up a informational advisor to contain the current time-step for
  // reference.  Activating this controller can also be used to freeze the
  // time-step at the current value.

  std::shared_ptr<fixed_ts_advisor> sp_dt(
      new fixed_ts_advisor("Current Time-Step", ts_advisor::req, dt, false));
  mngr.add_advisor(sp_dt);

  // Set up a required time-step to be activated at the user's discretion

  std::shared_ptr<fixed_ts_advisor> sp_ovr(new fixed_ts_advisor(
      "User Override", ts_advisor::req, override_dt, false));
  mngr.add_advisor(sp_ovr);

  // Set up a min timestep

  std::shared_ptr<fixed_ts_advisor> sp_min(
      new fixed_ts_advisor("Minimum", ts_advisor::min, ts_advisor::ts_small()));
  mngr.add_advisor(sp_min);
  sp_min->set_fixed_value(dt_min);

  // Set up a lower limit on the timestep rate of change

  std::shared_ptr<ratio_ts_advisor> sp_llr(
      new ratio_ts_advisor("Rate of Change Lower Limit", ts_advisor::min, 0.8));
  mngr.add_advisor(sp_llr);

  // Test the accessor
  {
    Remember(double tmp(sp_llr->get_ratio()););
    Check(rtt_dsxx::soft_equiv(tmp, 0.8));
  }

  // Set up an upper limit on the time-step rate of change

  std::shared_ptr<ratio_ts_advisor> sp_ulr(
      new ratio_ts_advisor("Rate of Change Upper Limit"));
  mngr.add_advisor(sp_ulr);

  // Set up an advisor to watch for an upper limit on the time-step.

  std::shared_ptr<fixed_ts_advisor> sp_max(new fixed_ts_advisor("Maximum"));
  mngr.add_advisor(sp_max);
  sp_max->set_fixed_value(dt_max);

  // Set up a target time advisor

  std::shared_ptr<target_ts_advisor> sp_gd(
      new target_ts_advisor("Graphics Dump", ts_advisor::max, graphics_time));
  mngr.add_advisor(sp_gd);

  // Now that all the advisors have been set up, perform time cycles

  for (int i = icycle_first; i <= icycle_last; ++i) {

    time = time + dt; //end of cycle time
    mngr.set_cycle_data(dt, i, time);

    // Make any user directed changes to controllers

    sp_dt->set_fixed_value(dt);
    if (override_flag) {
      // This branch is never used.  The next functions are commented out.
      // If you get here, please add a unit test and re-enable these
      // functions.
      ITFAILS;
      // sp_ovr->activate();
      // sp_ovr->set_fixed_value(override_dt);
    } else {
      sp_ovr->deactivate();
    }

    // Pass in the advisors owned by package_XXX for
    // that package to update

    xxx.advance_state();

    // Compute a new time-step and print results to screen

    dt = mngr.compute_new_timestep();
    mngr.print_summary();
  }

  // Dump a list of the advisors to the screen

  mngr.print_advisors();

  // Dump the advisor states for visual examination.

  mngr.print_adv_states();

  //------------------------------------------------------------//
  // Confirm that at least some of the output is correct.
  //------------------------------------------------------------//

  // Reference Values:
  double const prec(1.0e-5);
  double const ref1(3.345679);
  double const ref2(1.234568);
  double const ref3(1.371742);
  double const ref4(1.000000e-06);
  double const ref5(9.876543e-01);
  double const ref6(1.0);
  double const ref7(1.234568);
  double const ref8(1.481481);
  double const ref9(6.654321);
  double const ref10(1.000000e+05);
  double const ref11(1.371742);
  double const ref12(1.496914);
  double const ref13(2.716049);

  // Check final values:
  // ------------------------------
  if (mngr.get_cycle() == icycle_last)
    PASSMSG("get_cycle() returned the expected cycle index.");
  else
    FAILMSG("get_cycle() failed to return the expected cycle index.");

  if (mngr.get_controlling_advisor() == "Electron Temperature")
    PASSMSG("get_controlling_advisor() returned expected string.");
  else
    FAILMSG("get_controlling_advisor() failed to return the expected string.");

  if (rtt_dsxx::soft_equiv(ref11, xxx.get_dt_rec_te(), prec))
    PASSMSG("get_dt_rec_te() gave expected value.");
  else
    FAILMSG("get_dt_rec_te() did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref12, xxx.get_dt_rec_ti(), prec))
    PASSMSG("get_dt_rec_ti() gave expected value.");
  else
    FAILMSG("get_dt_rec_ti() did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref13, xxx.get_dt_rec_ri(), prec))
    PASSMSG("get_dt_rec_ri() gave expected value.");
  else
    FAILMSG("get_dt_rec_ri() did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref1, mngr.get_time(), prec))
    PASSMSG("get_time() gave expected value.");
  else
    FAILMSG("get_time() did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref2, mngr.get_dt(), prec))
    PASSMSG("get_dt() gave expected value.");
  else
    FAILMSG("get_dt() did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref3, mngr.get_dt_new(), prec))
    PASSMSG("get_dt_new() gave expected value.");
  else
    FAILMSG("get_dt_new() did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref4, sp_min->get_dt_rec(mngr), prec))
    PASSMSG(" sp_min->get_dt_rec(mngr) gave expected value.");
  else
    FAILMSG("sp_min->get_dt_rec(mngr) did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref5, sp_llr->get_dt_rec(mngr), prec))
    PASSMSG(" sp_llr->get_dt_rec(mngr) gave expected value.");
  else
    FAILMSG("sp_llr->get_dt_rec(mngr) did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref6, sp_ovr->get_dt_rec(mngr), prec))
    PASSMSG(" sp_ovr->get_dt_rec(mngr) gave expected value.");
  else
    FAILMSG("sp_ovr->get_dt_rec(mngr) did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref7, sp_dt->get_dt_rec(mngr), prec))
    PASSMSG(" sp_dt->get_dt_rec(mngr) gave expected value.");
  else
    FAILMSG("sp_dt->get_dt_rec(mngr) did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref8, sp_ulr->get_dt_rec(mngr), prec))
    PASSMSG(" sp_ulr->get_dt_rec(mngr) gave expected value.");
  else
    FAILMSG("sp_ulr->get_dt_rec(mngr) did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref9, sp_gd->get_dt_rec(mngr), prec))
    PASSMSG(" sp_gd->get_dt_rec(mngr) gave expected value.");
  else
    FAILMSG("sp_gd->get_dt_rec(mngr) did not give expected value.");

  if (rtt_dsxx::soft_equiv(ref10, sp_max->get_dt_rec(mngr), prec))
    PASSMSG(" sp_max->get_dt_rec(mngr) gave expected value.");
  else
    FAILMSG("sp_max->get_dt_rec(mngr) did not give expected value.");

  // Test the set_ratio manipulator for the ratio_ts_advisor
  // - reset the ratio_value to the default (1.2).
  {
    double const default_ratio(1.2);
    sp_llr->set_ratio();
    if (rtt_dsxx::soft_equiv(sp_llr->get_ratio(), default_ratio))
      PASSMSG("ratio_ts_advisor set_ratio manipultor/accessors work.");
    else
      FAILMSG("ratio_ts_advisor set_ratio manipultor/accessors are failing.");
  }

  // Test the set_target manipulator for the target_ts_advisor
  {
    double const ref_val(100.0);
    sp_gd->set_target(ref_val);
    if (rtt_dsxx::soft_equiv(sp_gd->get_target(), ref_val))
      PASSMSG("target_ts_advisor set_target manipultor/accessors work.");
    else
      FAILMSG("target_ts_advisor set_target manipultor/accessors are failing.");
  }

  // Check to make sure all processes passed.

  int npassed = (ut.numPasses > 0 && ut.numFails == 0) ? 1 : 0;
  rtt_c4::global_sum(npassed);

  if (npassed == rtt_c4::nodes()) {
    PASSMSG("All tests passed on all procs.");
  } else {
    std::ostringstream msg;
    msg << "Some tests failed on processor " << rtt_c4::node() << std::endl;
    FAILMSG(msg.str());
  }

  return;
}

//---------------------------------------------------------------------------//

void check_field_ts_advisor(rtt_dsxx::UnitTest &ut) {
  std::cout << "\nChecking the field_ts_advisor class...\n" << std::endl;

  rtt_timestep::field_ts_advisor ftsa;

  // Check manipulators
  std::cout << "Setting Frac Change to 1.0..." << std::endl;
  ftsa.set_fc(1.0);
  std::cout << "Setting Floor Value to 0.0001..." << std::endl;
  ftsa.set_floor(0.0001);
  std::cout << "Setting Update Method to q_mean..." << std::endl;
  ftsa.set_update_method(rtt_timestep::field_ts_advisor::q_mean);

  // Dump the state to an internal buffer and inspect the results
  std::ostringstream msg;
  ftsa.print_state(msg);
  std::cout << msg.str() << std::endl;

  { // Check the Fraction Change value
    std::string const expected("Fract Change   : 1");

    // find the line of interest
    std::string output(msg.str());
    size_t beg(output.find("Fract Change"));
    if (beg == std::string::npos) {
      FAILMSG("Did not find expected string!");
      return;
    }
    size_t end(output.find_first_of("\n", beg));
    if (beg == std::string::npos) {
      FAILMSG("Did not find expected string!");
      return;
    }
    std::string line(output.substr(beg, end - beg));
    if (line == expected)
      PASSMSG("'Fract Change' was set correctly.");
    else
      FAILMSG("Failed to set 'Fract Change' correctly.");
  }

  { // Check the Floor value
    std::string const expected("Floor Value    : 0.0001");

    // find the line of interest
    std::string output(msg.str());
    size_t beg(output.find("Floor Value"));
    if (beg == std::string::npos) {
      FAILMSG("Did not find expected string!");
      return;
    }
    size_t end(output.find_first_of("\n", beg));
    if (beg == std::string::npos) {
      FAILMSG("Did not find expected string!");
      return;
    }
    std::string line(output.substr(beg, end - beg));
    if (line == expected)
      PASSMSG("'Floor Value' was set correctly.");
    else
      FAILMSG("Failed to set 'Floor Value' correctly.");
  }

  { // Check the Update Method value
    std::string const expected("Update Method  : weighted by field value");

    // find the line of interest
    std::string output(msg.str());
    size_t beg(output.find("Update Method"));
    if (beg == std::string::npos) {
      FAILMSG("Did not find expected string!");
      return;
    }
    size_t end(output.find_first_of("\n", beg));
    if (beg == std::string::npos) {
      FAILMSG("Did not find expected string!");
      return;
    }
    std::string line(output.substr(beg, end - beg));
    if (line == expected)
      PASSMSG("'Update Method' was set correctly.");
    else
      FAILMSG("Failed to set 'Update Method' correctly.");
  }

  return;
}

//---------------------------------------------------------------------------//
// end of tstTimeStep.c
//---------------------------------------------------------------------------//
