//----------------------------------*-C++-*----------------------------------//
// target_ts_advisor.cc
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> Defines the target time-step advisor.
//---------------------------------------------------------------------------//

#include "timestep/target_ts_advisor.hh"

#include "timestep/ts_manager.hh"

#include "ds++/Assert.hh"

#include "c4/global.hh"

#include <iostream>

using std::cout;
using std::endl;

using namespace rtt_timestep;

target_ts_advisor::target_ts_advisor(
    const std::string   &name_,
    const usage_flag usage_,
    const double target_value_,
    const bool active_)

    : ts_advisor(name_, usage_, active_),
      target_value(target_value_)

{
    Ensure(invariant_satisfied());
}

target_ts_advisor::~target_ts_advisor()
{
// empty
}

double target_ts_advisor::get_dt_rec(const ts_manager &tsm) const
{
    Require(invariant_satisfied());

    double dt_rec = target_value - tsm.get_time();
    if (dt_rec <= small())
    {
	dt_rec = large();
    }
    return dt_rec;
}

void target_ts_advisor::print_state() const
{
    if (C4::node() != 0)
	return;
    
    std::string status = is_active() ? "true " : "false";
    cout << endl;
    cout << "  ** Time-Step Advisor State Listing **" << endl;
    cout << "  Name - " << get_name() << endl;
    cout << "  Type           : " << "Target Advisor" << endl;
    cout << "  Active         : " << status << endl;
    cout << "  Usage          : " << usage_flag_name(get_usage()) << endl;
    cout << "  Target Value   : " << target_value << endl;
    cout << endl;
}

bool target_ts_advisor::invariant_satisfied() const
{
    return ts_advisor::invariant_satisfied();
}

//---------------------------------------------------------------------------//
//                              end of target_ts_advisor.cc
//---------------------------------------------------------------------------//
