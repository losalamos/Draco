//----------------------------------*-C++-*----------------------------------//
// ratio_ts_advisor.cc
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> Defines the ratio time-step advisor.
//---------------------------------------------------------------------------//

#include "timestep/ratio_ts_advisor.hh"

#include "ds++/Assert.hh"

#include <iostream>

using std::cout;
using std::endl;

ratio_ts_advisor::ratio_ts_advisor(const std::string &name_,
				   const usage_flag usage_,
				   const double ratio_value_,
				   const bool active_) 

    : ts_advisor(name_, usage_, active_), ratio_value(ratio_value_)

{
    Ensure(invariant_satisfied());
}


ratio_ts_advisor::~ratio_ts_advisor()
{
// empty
}

void ratio_ts_advisor::update_tstep(const double current_dt,
				    const int cycle_)

{
    Require(invariant_satisfied());
    Require(current_dt > 0.)
	dt_rec = ratio_value*current_dt;
    cycle_at_last_update = cycle_;
    Ensure(invariant_satisfied());
}

void ratio_ts_advisor::print_state() const
{
    std::string status = active ? "true " : "false";
    cout << endl;
    cout << "  ** Time-Step Advisor State Listing **" << endl;
    cout << "  Name - " << name << endl;
    cout << "  Type           : " << "Ratio Advisor" << endl;
    cout << "  Active         : " << status << endl;
    cout << "  Usage          : " << usage_flag_name(usage) << endl;
    cout << "  Last Update    : " << "cycle " << cycle_at_last_update << endl;
    cout << "  Ratio Value    : " << ratio_value << endl;
    cout << "  dt_rec         : " << dt_rec << endl;
    cout << endl;
}

bool ratio_ts_advisor::invariant_satisfied() const
{
    bool ldum =
	name.length() != 0 &&
	inf    <= usage &&
	usage  <= req &&
	0. < dt_rec &&
        0. < ratio_value;

    return ldum;
}

//---------------------------------------------------------------------------//
//                              end of ratio_ts_advisor.cc
//---------------------------------------------------------------------------//
