//----------------------------------*-C++-*----------------------------------//
// field_ts_advisor.hh
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> Defines the field time-step advisor.
//---------------------------------------------------------------------------//

#ifndef __timestep_field_ts_advisor_hh__
#define __timestep_field_ts_advisor_hh__

//===========================================================================//
// class field_ts_advisor - Estimates a new timestep based on current fields.
//
// This class provides a means to estimate a
// new timestep based on the current dQ/dt where Q is some field quantity
// to be controlled (i.e. temperature, energy, particle number density,
// etc....), and dt is the current timestep.  
// 
//===========================================================================//

#include "timestep/ts_advisor.hh"
#include <iostream>

using std::cerr;

namespace rtt_timestep {

class field_ts_advisor : public ts_advisor {

// NESTED CLASSES AND TYPEDEFS

  public:

// Flag to determine the method used to produce the recommended 
// timestep. The recommended timestep will be:
// Based on a norm of a control function "alpha", where
// alpha = abs(del_Q/Q_old), where Q is a field  being monitored,
// i.e. temperature, energy, particles, etc. 
// Alpha is computed point-wise in the field, then a vector
// norm is applied to the point-wise values. The available 
// norms are:
 
    enum update_method_flag {

	inf_norm,   // infinity norm (max)
	a_mean,     // arithmetic mean
	q_mean,     // Q weighted mean
	rc_mean,    // relative change (alpha) weighted mean
	rcq_mean,   // product of Q and alpha weighted mean
        last_umf    // dummy to mark end of list
    };    

// DATA

  private:

    update_method_flag update_method; //update method for dt_rec
    double fc_value;                  //frac change  value for field advisor
    double floor_value;               //floor value for field advisor
    int    cycle_at_last_update;      //problem time-cycle index at last update
    double dt_rec;                    //the recommended time-step

// STATIC CLASS METHODS

  public:

    static std::string update_method_flag_name(const int i)
    {
	static const std::string update_method_flag_names [last_umf] =
	{ "infinity norm",
	  "arithmetic mean",
	  "weighted by field value",
	  "weighted by relative change",
	  "field value and relative change"
	};
	return update_method_flag_names[i];
    };

// CREATORS

    field_ts_advisor(const std::string &name_ = std::string("Unlabeled"),
		     const usage_flag usage_ = max,
		     const update_method_flag update_method_ = inf_norm,
		     const double fc_value_ = 0.1,
		     const double floor_value_ = small(),
		     const bool active_ = true);
    
    ~field_ts_advisor();

// MANIPULATORS

// A utility function to calculate a floor as a
// fraction of the max value in a field

    template < class FT >
    void set_floor(const FT &y1, double frac=0.001);

// Update the recommended time-step for advisors based on fields.
// q_old is the field value at the beginning of the current time-step,
// q_new is the field value at the end of the current time-step.

    template < class FT >
    void update_tstep(const ts_manager &tsm,
		      const FT &q_old, 
		      const FT &q_new);

// Set the fractional change value

    void set_fc(const double value_ = 0.1)
    { 
	fc_value = value_;
    }

// Set the floor value 

    void set_floor(const double value_ = small() )
    { 
	floor_value = value_;
    }

// Set update method

    void set_update_method(const update_method_flag  flag = inf_norm)
    {
	update_method = flag;
    }

// ACCESSORS

// Produce the recommended time-step

    double get_dt_rec(const ts_manager &tsm) const;

// Determine if the advisor is fit to use in
// a time-step calculation

    bool advisor_usable(const ts_manager &tsm) const;

// Print state

    void print_state() const;

// Invariant function

    bool invariant_satisfied() const;

};

} //end namespace rtt_timestep

#endif                          // __timestep_field_ts_advisor_hh__

//---------------------------------------------------------------------------//
//                              end of timestep/field_ts_advisor.hh
//---------------------------------------------------------------------------//
