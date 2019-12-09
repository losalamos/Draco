//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   dummy_package.cc
 * \author John McGhee
 * \date   Thu Aug 27 07:48:41 1998
 * \brief  A dummy package to exercize the time-step controller field
 *         advisors.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#include "dummy_package.hh"
#include "../field_ts_advisor.hh"
#include "../ts_manager.hh"
#include <vector>

namespace rtt_timestep_test {

dummy_package::dummy_package(rtt_timestep::ts_manager &tsm_)
    : tsm(tsm_), sp_te(new rtt_timestep::field_ts_advisor(
                     "Electron Temperature", rtt_timestep::ts_advisor::max,
                     rtt_timestep::field_ts_advisor::a_mean)),
      sp_ti(new rtt_timestep::field_ts_advisor(
          "Ion Temperature", rtt_timestep::ts_advisor::max,
          rtt_timestep::field_ts_advisor::rc_mean)),
      sp_ri(new rtt_timestep::field_ts_advisor(
          "Radiation Intensity", rtt_timestep::ts_advisor::max,
          rtt_timestep::field_ts_advisor::rcq_mean)) {
  // Set up a Electron Temperature advisor
  tsm.add_advisor(sp_te);

  // Set up a Ion-Temperature advisor
  tsm.add_advisor(sp_ti);

  // Set up a Radiation-Intensity advisor
  tsm.add_advisor(sp_ri);
}

dummy_package::~dummy_package() {
  tsm.remove_advisor(sp_te);
  tsm.remove_advisor(sp_ti);
  tsm.remove_advisor(sp_ri);
}

void dummy_package::advance_state() {
  /*!
     * \bug RMR The following code does not adequately cover all of the
     *      possible failure modes.  Please update this method with fuller
     *      coverage. 
     */

  // Create a set of dummy arrays to serve as control fields for use in
  // exercizing the various advisors.

  const double a1Array[] = {1., 10., 11., 3., 2., 5., 5., 6.7};
  const int sizeaArray = sizeof(a1Array) / sizeof(a1Array[0]);

  const int sizea = sizeaArray;
  const double *a1 = a1Array;

  std::vector<double> te_old(a1, a1 + sizea);
  std::vector<double> te_new = element_wise_multiply(1.09, te_old);
  std::vector<double> ti_old = element_wise_multiply(0.97, te_old);
  std::vector<double> ti_new = element_wise_multiply(1.05, te_old);
  std::vector<double> ri_old = element_wise_multiply(1.10, te_old);
  std::vector<double> ri_new = element_wise_multiply(1.15, te_old);

  // Set a floor for the electron temperature controller, to execcize this
  // method. Just accelpt the default floor on the other controllers.

  sp_te->set_floor(te_new, 0.001);

  // Get a new time-step from each of the advisors that belong to this
  // package.

  sp_te->update_tstep(tsm, te_old, te_new);
  sp_ti->update_tstep(tsm, ti_old, ti_new);
  sp_ri->update_tstep(tsm, ri_old, ri_new);

  return;
}

//---------------------------------------------------------------------------//

std::vector<double>
dummy_package::element_wise_multiply(double const a,
                                     std::vector<double> const &v) {
  std::vector<double> results(v.size());

  for (size_t i = 0; i < v.size(); ++i)
    results[i] = a * v[i];

  return results;
}

//---------------------------------------------------------------------------//

//! \brief Provide recommended dt values (TE).
double dummy_package::get_dt_rec_te() const { return sp_te->get_dt_rec(tsm); }

//! \brief Provide recommended dt values (TI).
double dummy_package::get_dt_rec_ti() const { return sp_ti->get_dt_rec(tsm); }

//! \brief Provide recommended dt values (RI).
double dummy_package::get_dt_rec_ri() const { return sp_ri->get_dt_rec(tsm); }

} // end namespace rtt_timestep_test

//---------------------------------------------------------------------------//
// end of dummy_package.cc
//---------------------------------------------------------------------------//
