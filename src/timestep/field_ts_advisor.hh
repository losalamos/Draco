//----------------------------------*-C++-*----------------------------------//
// field_ts_advisor.hh
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __timestep_field_ts_advisor_hh__
#define __timestep_field_ts_advisor_hh__

//===========================================================================//
// class field_ts_advisor - This class provides a means to estimate a
// new timestep based on the current dQ/dt where Q is some field quantity
// to be controlled (i.e. temperature, energy, particle number density,
// etc....), and dt is the current timestep.  
// 
//===========================================================================//

#include "timestep/ts_advisor.hh"

class field_ts_advisor : public ts_advisor {
    
  public:

// Flag to determine the method used to produce the recommended 
// timestep. The recommended timestep will be:
// Based on a norm of a control function "alpha", where
// alpha = abs(del_Q/Q), where Q is a field  being monitored,
// i.e. temperature, energy, particles, etc. The available 
// norms are:
 
    enum update_method_flag {

	inf_norm,   // infinity norm (max)
	a_mean,     // arithmetic mean
	q_mean,     // Q weighted mean
	rc_mean,    // relative change (alpha) weighted mean
	rcq_mean    // product of Q and alpha weighted mean
    };

    static std::string update_method_flag_name(const int i)
    {
	static const std::string update_method_flag_names [5] =
	{ "infinity norm",
	  "arithmetic mean",
	  "weighted by field value",
	  "weighted by relative change",
	  "field value and relative change"
	};
	return update_method_flag_names[i];
    };

    field_ts_advisor(const std::string &name_ = std::string("Unlabeled"),
		     const usage_flag usage_ = max,
		     const update_method_flag update_method_ = inf_norm,
		     const double fc_value_ = 0.1,
		     const double floor_value_ = small(),
		     const bool active_ = true);
    
    ~field_ts_advisor();

// Update the recommended time-step for advisors based on fields

    template < class FT >
    void update_tstep(const FT &y1, const FT &y2, 
		      double current_dt, 
		      int cycle_);

// A utility function to calculate a floor as a
// fraction of the max value in a field

    template < class FT >
    void set_floor(const FT &y1, double frac=0.001);

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

// Print state

    void print_state();

  private:

    bool invariant_satisfied();       //invariant function
    update_method_flag update_method; //update method for dt_rec
    double fc_value;                  //frac change  value for field advisor
    double floor_value;               //floor value for field advisor

};

#endif                          // __timestep_field_ts_advisor_hh__

//---------------------------------------------------------------------------//
//                              end of timestep/field_ts_advisor.hh
//---------------------------------------------------------------------------//
