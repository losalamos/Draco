//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ts_advisor.cc
 * \author John McGhee
 * \date   Thu Apr  2 14:06:18 1998
 * \brief  Defines the base class time-step advisor.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ts_advisor.hh"
#include "c4/C4_Functions.hh"
#include <string>

namespace rtt_timestep {

//----------------------------------------------------------------------------//
ts_advisor::ts_advisor(const std::string &name_, const usage_flag usage_,
                       const bool active_)
    : name(name_), usage(usage_), active(active_) {
  // empty
}

//----------------------------------------------------------------------------//
void ts_advisor::print(const ts_manager &tsm, const bool controlling) const {
  using std::cout;
  using std::endl;
  using std::string;

  if (rtt_c4::node() != 0)
    return;

  string status = advisor_usable(tsm) ? "true " : "false";
  string space = "   ";
  string cflag = controlling ? "  ==> " : "      ";
  cout << cflag;
  cout.width(12);
  cout << get_dt_rec(tsm) << space << status << space << name << endl;
}

} // namespace rtt_timestep

//---------------------------------------------------------------------------//
// end of ts_advisor.cc
//---------------------------------------------------------------------------//
