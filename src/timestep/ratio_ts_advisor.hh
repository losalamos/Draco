//----------------------------------*-C++-*----------------------------------//
// ratio_ts_advisor.hh
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> Defines the ratio time-step advisor.
//---------------------------------------------------------------------------//

#ifndef __timestep_ratio_ts_advisor_hh__
#define __timestep_ratio_ts_advisor_hh__

//===========================================================================//
// class ratio_ts_advisor - Calculates a new timestep as a ratio of 
//                          the current time-step
//
// This class provides a means to calculate a
// new timestep as a ratio of the current time-step. This is useful to
// limit the rate of change in the time-step from one time cycle to the 
// next.
// 
//===========================================================================//

#include "ts_advisor.hh"

namespace rtt_timestep {

class ratio_ts_advisor : public ts_advisor {


// DATA

 // Value used to compute ratio advisor

  private:
    double ratio_value; 



// CREATORS
  public:

    ratio_ts_advisor(const std::string &name_  = std::string("Unlabeled"),
		     const usage_flag usage_ = max,
		     const double ratio_value_ = 1.20,
		     const bool active_ = true);
    
    ~ratio_ts_advisor();



// MANIPULATORS


// Set the ratio value

    void set_ratio(const double value_ = 1.2)
    { 
	ratio_value = value_;
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

#endif                          // __timestep_ratio_ts_advisor_hh__

//---------------------------------------------------------------------------//
//                              end of ratio_ts_advisor.hh
//---------------------------------------------------------------------------//
