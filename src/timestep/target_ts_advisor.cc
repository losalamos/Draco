//----------------------------------*-C++-*----------------------------------//
// target_ts_advisor.cc
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "timestep/target_ts_advisor.hh"

#include "ds++/Assert.hh"

#include <iostream>

using std::cout;
using std::endl;

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

void target_ts_advisor::update_tstep( const double end_of_cycle_time,
				      const int cycle_)
{
    Require(invariant_satisfied());
    
    dt_rec = target_value - end_of_cycle_time;
    if (dt_rec <= small ())
    {
	dt_rec = large();
    }
    
    cycle_at_last_update = cycle_;
    Ensure(invariant_satisfied());
}


bool target_ts_advisor::invariant_satisfied()
{
    bool ldum =
	name.length() != 0 &&
	inf    <= usage &&
	usage  <= req &&
	0. < dt_rec;

    return ldum;
}

void target_ts_advisor::print_state()
{
    std::string status = active ? "true " : "false";
    cout << endl;
    cout << "  ** Time-Step Advisor State Listing **" << endl;
    cout << "  Name - " << name << endl;
    cout << "  Type           : " << "Target Advisor" << endl;
    cout << "  Active         : " << status << endl;
    cout << "  Usage          : " << usage_flag_name(usage) << endl;
    cout << "  Last Update    : " << "cycle " << cycle_at_last_update << endl;
    cout << "  Target Value   : " << target_value << endl;
    cout << "  dt_rec         : " << dt_rec << endl;
    cout << endl;
}
//---------------------------------------------------------------------------//
//                              end of target_ts_advisor.cc
//---------------------------------------------------------------------------//
