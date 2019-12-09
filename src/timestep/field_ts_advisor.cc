//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   field_ts_advisor.cc
 * \author John McGhee
 * \date   Thu Apr  2 14:06:18 1998
 * \brief  Defines the field time-step advisor.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.  
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "field_ts_advisor.hh"
#include "ts_manager.hh"
#include "c4/C4_Functions.hh"

using std::cout;
using std::endl;

namespace rtt_timestep {

field_ts_advisor::field_ts_advisor(const std::string &name_,
                                   const usage_flag usage_,
                                   const update_method_flag update_method_,
                                   const double fc_value_,
                                   const double floor_value_,
                                   const bool active_)

    : ts_advisor(name_, usage_, active_), update_method(update_method_),
      fc_value(fc_value_), floor_value(floor_value_),
      cycle_at_last_update(-989898), dt_rec(large()) {
  Ensure(invariant_satisfied());
}

// Produce the recommended time-step

double field_ts_advisor::get_dt_rec(const ts_manager &tsm) const {
  using std::cerr;

  if (cycle_at_last_update != tsm.get_cycle() && is_active() &&
      rtt_c4::node() == 0) {
    cerr << "Warning: ts_adivsor " << get_name() << " has not been updated"
         << endl;
    cerr << "         and will not be used in time-step calculations." << endl;
  }
  return dt_rec;
}

bool field_ts_advisor::advisor_usable(const ts_manager &tsm) const {
  return (is_active()) && (cycle_at_last_update == tsm.get_cycle());
}

void field_ts_advisor::print_state(std::ostream &out) const {
  if (rtt_c4::node() != 0)
    return;

  std::string status = is_active() ? "true " : "false";
  out << endl;
  out << "  ** Time-Step Advisor State Listing **" << endl;
  out << "  Name - " << get_name() << endl;
  out << "  Type           : "
      << "Field Advisor" << endl;
  out << "  Active         : " << status << endl;
  out << "  Usage          : " << usage_flag_name(get_usage()) << endl;
  out << "  Last Update    : "
      << "cycle " << cycle_at_last_update << endl;
  out << "  Update Method  : " << update_method_flag_name(update_method)
      << endl;
  out << "  Fract Change   : " << fc_value << endl;
  out << "  Floor Value    : " << floor_value << endl;
  out << "  dt_rec         : " << dt_rec << endl;
  out << endl;
}

bool field_ts_advisor::invariant_satisfied() const {
  return ts_advisor::invariant_satisfied() && 0. < dt_rec && 0. < floor_value &&
         0. < fc_value && 0 <= update_method && update_method < last_umf;
}

} // namespace rtt_timestep

//---------------------------------------------------------------------------//
// end of field_ts_advisor.cc
//---------------------------------------------------------------------------//
