//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   timestep/fixed_ts_advisor.cc
 * \author John McGhee
 * \date   Mon Apr  6 17:22:53 1998
 * \brief  Defines the fixed time-step advisor.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "fixed_ts_advisor.hh"
#include "c4/C4_Functions.hh"
#include "ds++/Assert.hh"

namespace rtt_timestep {

//---------------------------------------------------------------------------//
fixed_ts_advisor::fixed_ts_advisor(const std::string &name_,
                                   const usage_flag usage_,
                                   const double fixed_value_,
                                   const bool active_)

    : ts_advisor(name_, usage_, active_), fixed_value(fixed_value_) {
  Ensure(invariant_satisfied());
}

//---------------------------------------------------------------------------//
double fixed_ts_advisor::get_dt_rec(const ts_manager & /*tsm*/) const {
  Require(invariant_satisfied());
  return fixed_value;
}

//---------------------------------------------------------------------------//
void fixed_ts_advisor::print_state(std::ostream &out) const {
  using std::endl;

  if (rtt_c4::node() != 0)
    return;

  std::string status = is_active() ? "true " : "false";
  out << endl;
  out << "  ** Time-Step Advisor State Listing **" << endl;
  out << "  Name - " << get_name() << endl;
  out << "  Type           : "
      << "Fixed Advisor" << endl;
  out << "  Active         : " << status << endl;
  out << "  Usage          : " << usage_flag_name(get_usage()) << endl;
  out << "  Fixed Value    : " << fixed_value << endl;
  out << endl;
}

//---------------------------------------------------------------------------//
bool fixed_ts_advisor::invariant_satisfied() const {
  return ts_advisor::invariant_satisfied() && 0. < fixed_value;
}

} // namespace rtt_timestep

//---------------------------------------------------------------------------//
// end of fixed_ts_advisor.cc
//---------------------------------------------------------------------------//
