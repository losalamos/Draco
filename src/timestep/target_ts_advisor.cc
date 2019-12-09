//----------------------------------*-C++-*----------------------------------//
/*!
 * \file    target_ts_advisor.cc
 * \author  John McGhee
 * \date    Thu Apr  2 14:06:18 1998
 * \brief   Defines the target time-step advisor.
 * \note    Copyright (C) 2016-2019 Triad National Security, LLC.
 *          All rights reserved.
 * \version $Id$
 */
//---------------------------------------------------------------------------//

#include "target_ts_advisor.hh"
#include "ts_manager.hh"
#include "c4/C4_Functions.hh"

namespace rtt_timestep {

target_ts_advisor::target_ts_advisor(const std::string &name_,
                                     const usage_flag usage_,
                                     const double target_value_,
                                     const bool active_)

    : ts_advisor(name_, usage_, active_), target_value(target_value_)

{
  Ensure(invariant_satisfied());
}

double target_ts_advisor::get_dt_rec(const ts_manager &tsm) const {
  Require(invariant_satisfied());

  double dt_rec = target_value - tsm.get_time();
  if (dt_rec <= ts_small()) {
    dt_rec = large();
  }
  return dt_rec;
}

void target_ts_advisor::print_state(std::ostream &out) const {
  using std::endl;
  if (rtt_c4::node() != 0)
    return;

  std::string status = is_active() ? "true " : "false";
  out << endl;
  out << "  ** Time-Step Advisor State Listing **" << endl;
  out << "  Name - " << get_name() << endl;
  out << "  Type           : "
      << "Target Advisor" << endl;
  out << "  Active         : " << status << endl;
  out << "  Usage          : " << usage_flag_name(get_usage()) << endl;
  out << "  Target Value   : " << target_value << endl;
  out << endl;
}

bool target_ts_advisor::invariant_satisfied() const {
  return ts_advisor::invariant_satisfied();
}

} // namespace rtt_timestep

//---------------------------------------------------------------------------//
// end of target_ts_advisor.cc
//---------------------------------------------------------------------------//
