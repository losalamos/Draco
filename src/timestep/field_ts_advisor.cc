//----------------------------------*-C++-*----------------------------------//
// field_ts_advisor.cc
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> Defines the field time-step advisor.
//---------------------------------------------------------------------------//

#include "timestep/field_ts_advisor.hh"

#include "ds++/Assert.hh"

#include <iostream>

using std::cout;
using std::endl;

field_ts_advisor::field_ts_advisor(const std::string &name_,
				   const usage_flag usage_,
				   const update_method_flag update_method_,
				   const double fc_value_, 
				   const double floor_value_,
				   const bool active_)

    : ts_advisor (name_, usage_, active_), 
      update_method(update_method_),
      fc_value(fc_value_), floor_value(floor_value_)
{
    Ensure(invariant_satisfied());
}


field_ts_advisor::~field_ts_advisor()
{
// empty
}



void field_ts_advisor::print_state() const
{
    std::string status = active ? "true " : "false";
    cout << endl;
    cout << "  ** Time-Step Advisor State Listing **" << endl;
    cout << "  Name - " << name << endl;
    cout << "  Type           : " << "Field Advisor" << endl;
    cout << "  Active         : " << status << endl;
    cout << "  Usage          : " << usage_flag_name(usage) << endl;
    cout << "  Last Update    : " << "cycle " << cycle_at_last_update << endl;
    cout << "  Update Method  : " << 
	update_method_flag_name(update_method) << endl;
    cout << "  Fract Change   : " << fc_value << endl;
    cout << "  Floor Value    : " << floor_value << endl;
    cout << "  dt_rec         : " << dt_rec << endl;
    cout << endl;
}

bool field_ts_advisor::invariant_satisfied() const
{
    bool ldum =
	name.length() != 0 &&
	0      <= usage &&
	usage  <  last_usage  &&
	0. < dt_rec &&
	0. < floor_value  &&
	0. < fc_value  &&
	0 <= update_method && 
	update_method < last_umf;

    return ldum;
}

//---------------------------------------------------------------------------//
//                              end of field_ts_advisor.cc
//---------------------------------------------------------------------------//
