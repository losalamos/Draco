//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ts_manager.cc
 * \author John McGhee
 * \date   Mon Apr  6 17:22:53 1998
 * \brief  Defines a manager utility for time-step advisors.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ts_manager.hh"
#include "c4/C4_Functions.hh"

using std::cerr;
using std::cout;
using std::endl;
using std::ios;
using std::list;

namespace rtt_timestep {

ts_manager::ts_manager(void)
    : dt_new(ts_advisor::ts_small()), time(0.0), dt(ts_advisor::ts_small()),
      cycle(9999), controlling_advisor("Not Set")

{
  Ensure(invariant_satisfied());
}

void ts_manager::add_advisor(const std::shared_ptr<ts_advisor> &new_advisor) {
  bool not_in_use = true;
  for (list<std::shared_ptr<ts_advisor>>::iterator py = advisors.begin();
       py != advisors.end(); py++) {
    if ((**py).get_name() == (*new_advisor).get_name()) {
      not_in_use = false;
      break;
    }
  }
  if (not_in_use) {
    advisors.push_front(new_advisor);
  } else {
    throw std::runtime_error("Name for requested advisor already in use");
  }
}

void ts_manager::remove_advisor(
    const std::shared_ptr<ts_advisor> &advisor_to_remove) {
  for (list<std::shared_ptr<ts_advisor>>::iterator py = advisors.begin();
       py != advisors.end(); py++) {
    if ((**py).get_name() == (*advisor_to_remove).get_name()) {
      advisors.erase(py);
      return;
    }
  }
  throw std::runtime_error("Unable to find requested advisor");
}

void ts_manager::set_cycle_data(double dt_, int cycle_, double time_) {
  dt = dt_;
  cycle = cycle_;
  time = time_;
}

double ts_manager::compute_new_timestep() {
  // Check to be sure that there is at least one usable advisor

  bool found = false;
  for (list<std::shared_ptr<ts_advisor>>::iterator py = advisors.begin();
       py != advisors.end(); py++) {
    if ((**py).advisor_usable(*this)) {
      found = true;
      break;
    }
  }

  if (!found) {
    if (rtt_c4::node() == 0) {
      cerr << "  ** Time-Step Manager Warning **" << endl;
      cerr << "  No usable time-step advisors found," << endl;
      cerr << "  defaulting to current time-step" << endl;
    }
    dt_new = dt;
    controlling_advisor = "Current Time-Step";
    return dt_new;
  }

  // Check for one and only one mandatory advisor

  int i = 0;
  for (list<std::shared_ptr<ts_advisor>>::iterator py = advisors.begin();
       py != advisors.end(); py++) {
    if ((**py).advisor_usable(*this) && (**py).get_usage() == ts_advisor::req) {
      i++;
      dt_new = (**py).get_dt_rec(*this);
      controlling_advisor = (**py).get_name();
    }
  }
  if (i == 1) {
    return dt_new;
  } else if (i != 0 && rtt_c4::node() == 0) {
    cerr << "  ** Time-Step Manager Warning **" << endl;
    cerr << "  Cycle Number: " << cycle << endl;
    cerr << "  More than one mandatory advisor found," << endl;
    cerr << "  defaulting to last found" << endl;
  }

  // Loop over the advisors finding the one that controls

  list<std::shared_ptr<ts_advisor>>::iterator py1 = advisors.end();
  list<std::shared_ptr<ts_advisor>>::iterator py2 = advisors.end();
  double x1 = ts_advisor::ts_small();
  double x2 = ts_advisor::large();
  for (list<std::shared_ptr<ts_advisor>>::iterator py = advisors.begin();
       py != advisors.end(); py++) {
    if ((**py).advisor_usable(*this)) {
      if ((**py).get_usage() == ts_advisor::min) {
        if ((**py).get_dt_rec(*this) > x1) {
          x1 = (**py).get_dt_rec(*this);
          py1 = py;
        }
      } else if ((**py).get_usage() == ts_advisor::max) {
        if ((**py).get_dt_rec(*this) < x2) {
          x2 = (**py).get_dt_rec(*this);
          py2 = py;
        }
      }
    }
  }

  if (py1 == advisors.end() && py2 == advisors.end()) {
    if (rtt_c4::node() == 0) {
      cerr << "  ** Time-Step Manager Warning **" << endl;
      cerr << "  Cycle Number: " << cycle << endl;
      cerr << "  No usable time-step advisors found," << endl;
      cerr << "  defaulting to current time-step" << endl;
    }
    dt_new = dt;
    controlling_advisor = "Current Time-Step";
    return dt_new;
  } else if (py1 == advisors.end()) {
    dt_new = x2;
    controlling_advisor = (**py2).get_name();
  } else if (py2 == advisors.end()) {
    dt_new = x1;
    controlling_advisor = (**py1).get_name();
  } else {
    if (x1 > x2) {
      dt_new = x1;
      controlling_advisor = (**py1).get_name();
      if (rtt_c4::node() == 0) {
        cerr << "  ** Time-Step Manager Warning **" << endl;
        cerr << "  Cycle Number: " << cycle << endl;
        cerr << "  No window between min and max advisors," << endl;
        cerr << "  defaulting to min recommended dt" << endl;
      }
    } else {
      dt_new = x2;
      controlling_advisor = (**py2).get_name();
    }
  }

  return dt_new;
}

//---------------------------------------------------------------------------//
/*!
 * Defines a functor which determines if one timestep advisor is less than
 * another. This is done by comparing the recommended time step of each advisor.
 */
class sptsa_less_than {
public:
  // The timestep manager
  ts_manager const &tsm;

  // Constructs the sptsa_less_than functor
  sptsa_less_than(ts_manager const &tsm_) : tsm(tsm_) {
    // empty
  }

  // Defines the operator () for sptsa_less_than functor
  bool operator()(const std::shared_ptr<ts_advisor> &sp_lhs,
                  const std::shared_ptr<ts_advisor> &sp_rhs) const {
    return (sp_lhs->get_dt_rec(tsm) < sp_rhs->get_dt_rec(tsm));
  }
};

void ts_manager::print_advisors() const {
  if (rtt_c4::node() != 0)
    return;

  cout << endl;
  cout << "*** Time-Step Manager: Advisor Listing ***" << endl;
  for (list<std::shared_ptr<ts_advisor>>::const_iterator py = advisors.begin();
       py != advisors.end(); py++) {
    cout << (**py).get_name() << endl;
  }
  cout << endl;
}

void ts_manager::print_summary() const {
  if (rtt_c4::node() != 0)
    return;

  ios::fmtflags oldOptions = cout.flags();
  cout.setf(ios::scientific, ios::floatfield);
  cout.precision(4);
  cout.setf(ios::left, ios::adjustfield);
  list<std::shared_ptr<ts_advisor>> temp = advisors;
  temp.sort(sptsa_less_than(*this));
  cout << endl;
  cout << "  *** Time-Step Manager Summary ***" << endl;
  cout << "  Cycle Number         : " << cycle << endl;
  cout << "  Problem Time         : " << time << endl;
  cout << "  Current Time-Step    : " << dt << endl;
  cout << "  Recommended Time-Step: " << dt_new << endl;
  cout << "  Controlling Advisor  : " << controlling_advisor << endl;
  cout << endl;
  cout << "  *** Time-Step Advisor Data Table *** " << endl;
  cout << "    Recommendation   In-Use  Name " << endl;
  for (auto py = temp.begin(); py != temp.end(); py++) {
    (**py).print(*this, (**py).get_name() == controlling_advisor);
  }
  cout << endl;
  cout.flags(oldOptions);
}

void ts_manager::print_adv_states() const {
  if (rtt_c4::node() != 0)
    return;

  cout << endl;
  cout << "*** Time-Step Manager: Advisor State Listing ***" << endl;
  for (list<std::shared_ptr<ts_advisor>>::const_iterator py = advisors.begin();
       py != advisors.end(); py++) {
    (**py).print_state();
  }
  cout << endl;
}

bool ts_manager::invariant_satisfied() const {
  bool ldum = 0.0 < dt_new && 0.0 < dt && controlling_advisor.length() != 0;
  return ldum;
}

} // namespace rtt_timestep

//---------------------------------------------------------------------------//
// end of ts_manager.cc
//---------------------------------------------------------------------------//
