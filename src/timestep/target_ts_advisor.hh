//----------------------------------*-C++-*----------------------------------//
// target_ts_advisor.hh
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __timestep_target_ts_advisor_hh__
#define __timestep_target_ts_advisor_hh__

//===========================================================================//
// class target_ts_advisor - This class provides a means to calculate 
// a new time-step required to hit some "target" problem time. This
// is useful to assure that graphics dumps/IO/etc. occur precisely at
// a predetermined problem time.
// 
//===========================================================================//

#include "timestep/ts_advisor.hh"

class target_ts_advisor : public ts_advisor {


  public:
    
    target_ts_advisor( 
	const std::string &name_  = std::string("Unlabeled"),
	const usage_flag usage_ = max,
	const double target_value_ = -large(),
	const bool active_ = true );
    
    ~target_ts_advisor();

// Update the recommended time-step 
    
    void update_tstep( const double end_of_cycle_time,
		       const int cycle_ );

// Set the target value

    void set_target(const double value_ = -large() )
    { 
	target_value = value_;
    }

// Print state

    void print_state();

  private:

    bool invariant_satisfied();
    double target_value;              //target time (time)

};

#endif                          // __target_timestep_ts_advisor_hh__

//---------------------------------------------------------------------------//
//                              end of timestep/target_ts_advisor.hh
//---------------------------------------------------------------------------//
