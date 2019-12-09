//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   field_ts_advisor.t.hh
 * \author John McGhee
 * \date   Mon Aug 24 07:48:00 1998
 * \brief  Contains the template methods for the field ts_advisor class.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "field_ts_advisor.hh"
#include "ts_manager.hh"
#include "c4/C4_Functions.hh"
#include <cmath>

namespace rtt_timestep {

//---------------------------------------------------------------------------//
template <class FT>
void field_ts_advisor::set_floor(FT const &y1, double frac) {
  Require(invariant_satisfied());
  Require(frac > 0.);
  double x1 = -large();
  for (typename FT::const_iterator py1 = y1.begin(); py1 != y1.end(); py1++) {
    if (*py1 > x1) {
      x1 = *py1;
    }
  }
  x1 = x1 * frac;
  if (x1 <= ts_small()) {
    x1 = ts_small();
  }
  floor_value = x1;

  // All process will get the same floor_value.

  rtt_c4::global_min(floor_value);

  Ensure(invariant_satisfied());
}

//---------------------------------------------------------------------------//
template <class FT>
void field_ts_advisor::update_tstep(const ts_manager &tsm, const FT &q_old,
                                    const FT &q_new) {
  Require(invariant_satisfied());
  Require(tsm.get_dt() > 0.0);
  //    Require(FT::conformal(q_old,q_new));

  double x1 = 0.;
  double x2 = 0.;
  if (update_method == inf_norm) {
    x2 = 1.;
  }

  typename FT::const_iterator pq_new = q_new.begin();
  for (typename FT::const_iterator pq_old = q_old.begin();
       pq_old != q_old.end(); pq_old++, pq_new++) {

    if (*pq_old > floor_value && *pq_new > floor_value) {
      double delta_q = std::abs(*pq_new - *pq_old);
      double q_norm = *pq_old;
      double alpha = delta_q / q_norm;

      if (alpha < eps()) // Set noise to a hard zero
      {
        alpha = 0.;
        delta_q = 0.;
      }
      if (update_method == inf_norm) {
        if (alpha > x1)
          x1 = alpha;
      } else if (update_method == a_mean) {
        x1 = x1 + alpha;
        x2 = x2 + 1.;
      } else if (update_method == q_mean) {
        x1 = x1 + delta_q;
        x2 = x2 + q_norm;
      } else if (update_method == rc_mean) {
        x1 = x1 + alpha * alpha;
        x2 = x2 + alpha;
      } else if (update_method == rcq_mean) {
        x1 = x1 + alpha * delta_q;
        x2 = x2 + delta_q;
      } else {
        throw std::runtime_error("Unrecognized update method flag");
      }
    }
  }

  // If we are doing an infinity norm, then we must find the max on
  // all processors.
  // Otherwise we must add are partial sums across processors.

  if (update_method == inf_norm) {
    rtt_c4::global_max(x1);
  } else {
    rtt_c4::global_sum(x1);
    rtt_c4::global_sum(x2);
  }

  if (x1 < ts_small()) {
    dt_rec = large();
  } else {
    double fact = x2 * fc_value / x1;
    if (fact < ts_small()) {
      dt_rec = ts_small();
    } else {
      dt_rec = std::max(ts_small(), fact * tsm.get_dt());
    }
  }

  cycle_at_last_update = tsm.get_cycle();
  Ensure(invariant_satisfied());
}

} // namespace rtt_timestep

//---------------------------------------------------------------------------//
// end of field_ts_advisor.t.hh
//---------------------------------------------------------------------------//
