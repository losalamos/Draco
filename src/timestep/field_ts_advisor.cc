//----------------------------------*-C++-*----------------------------------//
// field_ts_advisor.cc
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "timestep/field_ts_advisor.hh"

#include "ds++/Assert.hh"

#include <algorithm>

#include <stdexcept>

#include <iostream>

#include <cmath>

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

template <class FT>
void field_ts_advisor::set_floor(const FT &y1, double frac)
{
    Require(invariant_satisfied());
    Require(frac > 0.);
    double x1 = -large();
    for (FT::const_iterator py1 = y1.begin(); py1 != y1.end(); py1++) 
    {
	if (*py1 > x1)
	{
	    x1 = *py1;
	}
    }
    x1 = x1*frac;
    if (x1 <= small())
    {
	x1 = small();
    }
    floor_value = x1;
    Ensure(invariant_satisfied());
}


template <class FT>
void field_ts_advisor::update_tstep(const FT &y1, const FT &y2, 
				    double current_dt,
				    int cycle_)
{

    Require(invariant_satisfied());
    Require(current_dt > 0.0);
//    Require(FT::conformal(y1,y2));

    double x1 = 0.;
    double x2 = 0.;
    if (update_method == inf_norm) 
    {
	x2 = 1.;
    }

    FT::const_iterator py2 = y2.begin();
    for (FT::const_iterator py1 = y1.begin(); py1 != y1.end(); py1++,py2++) 
    {
	
	if (*py1 > floor_value && *py2 > floor_value)
	{
	    double delta_y = std::abs(*py2-*py1);
	    double y_norm  = 0.5*(*py2+*py1);
	    double alpha = delta_y/y_norm;

	    if (alpha < eps()) // Set noise to a hard zero
	    {
		alpha = 0.;
		delta_y = 0.;
	    }
	    if ( update_method == inf_norm )
	    {
		if (alpha > x1) x1=alpha ;
	    }
	    else if ( update_method == a_mean )
	    {
		x1 = x1 + alpha;
		x2 = x2 + 1.;
	    }
	    else if ( update_method == q_mean )
	    {
		x1 = x1 + delta_y;
		x2 = x2 + y_norm;
	    }
	    else if ( update_method == rc_mean )
	    {
		x1 = x1 + alpha*alpha;
		x2 = x2 + alpha;
	    }
	    else if ( update_method == rcq_mean )
	    {
		x1 = x1 + alpha*delta_y;
		x2 = x2 + delta_y;
	    }
	    else
	    {
		throw std::runtime_error("Unrecognized update method flag");
	    }
	}
    }

    if (x1 < small()) 
    {
	dt_rec = large();
    }
    else 
    {
	double fact = x2*fc_value/x1;
	if (fact < small())
	{
	    dt_rec = small();
	}
	else
	{
	    dt_rec = std::max(small(),fact*current_dt);
	}
    }

    cycle_at_last_update = cycle_;
    Ensure(invariant_satisfied());

}

bool field_ts_advisor::invariant_satisfied()
{
    bool ldum =
	name.length() != 0 &&
	inf    <= usage &&
	usage  <= req &&
	0. < dt_rec &&
	0. < floor_value  &&
	0. < fc_value  &&
	inf_norm <= update_method && 
	update_method <= rcq_mean;

    return ldum;
}

void field_ts_advisor::print_state()
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



//---------------------------------------------------------------------------//
//                              end of field_ts_advisor.cc
//---------------------------------------------------------------------------//
