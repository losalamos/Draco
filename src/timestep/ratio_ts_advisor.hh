//----------------------------------*-C++-*----------------------------------//
// ratio_ts_advisor.hh
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __timestep_ratio_ts_advisor_hh__
#define __timestep_ratio_ts_advisor_hh__

//===========================================================================//
// class ratio_ts_advisor - This class provides a means to calculate a
// new timestep as a ratio of the current time-step. This is useful to
// limit the rate of change in the time-step form one time cycle to the 
// next.
// 
//===========================================================================//

#include "timestep/ts_advisor.hh"

class ratio_ts_advisor : public ts_advisor {

  public:

    ratio_ts_advisor(const std::string &name_  = std::string("Unlabeled"),
		     const usage_flag usage_ = max,
		     const double ratio_value_ = 1.20,
		     const bool active_ = true);
    
    ~ratio_ts_advisor();

// Update the recommended time-step 
    
    void update_tstep(const double current_dt,
		      const int cycle_);

// Set the ratio value

    void set_ratio(const double value_ = 1.2)
    { 
	ratio_value = value_;
    }

// Print state

    void print_state();

  private:

    bool invariant_satisfied();
    double ratio_value;  //value used to compute ratio advisor

};

#endif                          // __timestep_ratio_ts_advisor_hh__

//---------------------------------------------------------------------------//
//                              end of timestep/ratio_ts_advisor.hh
//---------------------------------------------------------------------------//
