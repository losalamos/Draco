//----------------------------------*-C++-*----------------------------------//
// fixed_ts_advisor.hh
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __timestep_fixed_ts_advisor_hh__
#define __timestep_fixed_ts_advisor_hh__

//===========================================================================//
// class fixed_ts_advisor - This class provides a means to introduce
// a user defined fixed value into the time-step calculation.
// This is useful to set min and max timesteps, or to force a
// timestep, etc.
// 
//===========================================================================//

#include "timestep/ts_advisor.hh"

class fixed_ts_advisor : public ts_advisor {

  public:

    fixed_ts_advisor( 
	const std::string &name_  = std::string("Unlabeled"),
	const usage_flag usage_ = max, 
	const double const_value_ = large(),
	const bool active_ = true);
    
    ~fixed_ts_advisor();

// Update the recommended time-step
    
    void update_tstep(const int cycle_);

// Set the fixed value

    void set_fixed_value(const double value_ = large())
    { 
	fixed_value = value_;
    }

// Print state

    void print_state();

  private:

    bool invariant_satisfied();

    double fixed_value;     //value used to oompute fixed advisor 

};

#endif                          // __timestep_fixed_ts_advisor_hh__

//---------------------------------------------------------------------------//
//                              end of timestep/fixed_ts_advisor.hh
//---------------------------------------------------------------------------//
