//----------------------------------*-C++-*----------------------------------//
// fixed_ts_advisor.cc
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> Defines the fixed time-step advisor. 
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "fixed_ts_advisor.hh"

#include "ds++/Assert.hh"

#include "c4/global.hh"

#include <iostream>

using std::cout;
using std::endl;

namespace rtt_timestep
{

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

double fixed_ts_advisor::get_dt_rec(const ts_manager &tsm) const
{
    Require(invariant_satisfied());
    return fixed_value;
}

void fixed_ts_advisor::print_state() const
{
    if (C4::node() != 0)
	return;
    
    std::string status = is_active() ? "true " : "false";
    cout << endl;
    cout << "  ** Time-Step Advisor State Listing **" << endl;
    cout << "  Name - " << get_name() << endl;
    cout << "  Type           : " << "Fixed Advisor" << endl;
    cout << "  Active         : " << status << endl;
    cout << "  Usage          : " << usage_flag_name(get_usage()) << endl;
    cout << "  Fixed Value    : " << fixed_value << endl;
    cout << endl;
}

bool fixed_ts_advisor::invariant_satisfied() const
{
    return ts_advisor::invariant_satisfied() &&
        0. < fixed_value;
}

} //end of rtt_timestep namespace

//---------------------------------------------------------------------------//
//                              end of fixed_ts_advisor.cc
//---------------------------------------------------------------------------//
