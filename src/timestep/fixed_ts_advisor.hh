//----------------------------------*-C++-*----------------------------------//
// fixed_ts_advisor.hh
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> Defines the fixed time-step advisor.
//---------------------------------------------------------------------------//

#ifndef __timestep_fixed_ts_advisor_hh__
#define __timestep_fixed_ts_advisor_hh__

//===========================================================================//
// class fixed_ts_advisor -  Introduces a user defined fixed value into the 
//                           time-step calculation.
//
// This class provides a means to introduce
// a user defined fixed value into the time-step calculation.
// This is useful to set min and max timesteps, or to force a
// timestep, etc.
// 
//===========================================================================//

#include "ts_advisor.hh"

namespace rtt_timestep {

class fixed_ts_advisor : public ts_advisor {

  // DATA

  // Value used to oompute fixed advisor 

  private:
    double fixed_value; 

   
// CREATORS

  public:

    fixed_ts_advisor( 
	const std::string &name_  = std::string("Unlabeled"),
	const usage_flag usage_ = max, 
	const double const_value_ = large(),
	const bool active_ = true);
    
    ~fixed_ts_advisor();

// MANIPULATORS


// Set the fixed value

    void set_fixed_value(const double value_ = large())
    { 
	fixed_value = value_;
    }


// ACCESSORS

// Update the recommended time-step
    
    double get_dt_rec(const ts_manager &tsm) const;

// Print state

    void print_state() const;

// Invariant function

    bool invariant_satisfied() const;

};

} // end of rtt_timestep namespace

#endif                          // __timestep_fixed_ts_advisor_hh__

//---------------------------------------------------------------------------//
//                              end of fixed_ts_advisor.hh
//---------------------------------------------------------------------------//
