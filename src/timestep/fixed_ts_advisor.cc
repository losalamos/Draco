//----------------------------------*-C++-*----------------------------------//
// fixed_ts_advisor.cc
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> Defines the fixed time-step advisor. 
//---------------------------------------------------------------------------//

#include "timestep/fixed_ts_advisor.hh"

#include "ds++/Assert.hh"

#include <iostream>

using std::cout;
using std::endl;

fixed_ts_advisor::fixed_ts_advisor(const std::string &name_,
				   const usage_flag usage_,
				   const double fixed_value_,
				   const bool active_)

    : ts_advisor(name_, usage_, active_),
      fixed_value(fixed_value_)
{
    Ensure(invariant_satisfied());
}


fixed_ts_advisor::~fixed_ts_advisor()
{
// empty
}

void fixed_ts_advisor::update_tstep(const int cycle_)
{
    Require(invariant_satisfied());
    dt_rec = fixed_value;
    cycle_at_last_update = cycle_;
    Ensure(invariant_satisfied());
}

void fixed_ts_advisor::print_state() const
{
    std::string status = active ? "true " : "false";
    cout << endl;
    cout << "  ** Time-Step Advisor State Listing **" << endl;
    cout << "  Name - " << name << endl;
    cout << "  Type           : " << "Fixed Advisor" << endl;
    cout << "  Active         : " << status << endl;
    cout << "  Usage          : " << usage_flag_name(usage) << endl;
    cout << "  Last Update    : " << "cycle " << cycle_at_last_update << endl;
    cout << "  Fixed Value    : " << fixed_value << endl;
    cout << "  dt_rec         : " << dt_rec << endl;
    cout << endl;
}

bool fixed_ts_advisor::invariant_satisfied() const
{
    bool ldum =
	name.length() != 0 &&
	inf    <= usage &&
	usage  <= req &&
	0. < dt_rec &&
        0. < fixed_value;

    return ldum;
}

//---------------------------------------------------------------------------//
//                              end of fixed_ts_advisor.cc
//---------------------------------------------------------------------------//
