//----------------------------------*-C++-*----------------------------------//
// ts_advisor.hh
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __timestep_ts_advisor_hh__
#define __timestep_ts_advisor_hh__

//===========================================================================//
// class ts_advisor - This is the base class time-step advisor.
// 
//===========================================================================//

#include <limits>

#include <string>

class ts_advisor {

  public:

    static double eps() 
    {
	return 100.*std::numeric_limits<double>::epsilon(); 
    }

    static double small() 
    {
	return 100.*std::numeric_limits<double>::min(); 
    }

    static double large() 
    {
	return 0.01*std::numeric_limits<double>::max(); 
    }

// Flag to determine how the recommended timestep is to be
// used. The recommended value "dt_rec" is to be considered:

    enum usage_flag {

        inf , // informational only, not to be used for control
	min , // a lower limit
	max , // a upper limit
	req   // a required value
    };

    static std::string usage_flag_name(const int i) 
    {
	static const std::string usage_flag_names [4] =
	{	"informational",
		"minimum",
		"maximum",
		"required"};
	return usage_flag_names[i];
    };

    ts_advisor(const std::string &name_  = std::string("Unlabeled"),
	       const usage_flag usage_ = max,
	       const bool active_ = true);
    
    virtual ~ts_advisor();

// Turn the advisor on and off

    void activate()
    {
	active = true;
    }

    void deactivate()
    {
	active = false;
    }

// Determine if the advisor is fit to use in
// a time-step calculation

    bool advisor_usable(const int cycle_) const
    {
	return (active == true) &&
	    (usage != inf) &&
	    (cycle_at_last_update == cycle_);
    }

// Get the usage

    usage_flag get_usage() const
    {
	return usage;
    }

// Get the name

    const std::string &get_name() const
    {
	return name;
    }
    
// Print out advisor data

    void print(const int cycle_, const bool controlling = false) const;

// Define the "less than" operator (<) so that the
// advisors can be sorted based on their recommended 
// time-steps

    bool operator<(const ts_advisor &rhs) const
    {
	return dt_rec < rhs.dt_rec;
    }

// Produce the recommended time-step

    double get_dt_rec()
    {
	return dt_rec;
    }

  protected:

    std::string name;                 //ID string
    usage_flag usage;                 //how to use dt_rec 
    bool   active;                    //on-off switch
    int    cycle_at_last_update;      //problem time-cycle index at last update
    double dt_rec;                    //the recommended time-step

};

#endif                          // __timestep_ts_advisor_hh__

//---------------------------------------------------------------------------//
//                              end of timestep/ts_advisor.hh
//---------------------------------------------------------------------------//
