//----------------------------------*-C++-*----------------------------------//
// ts_advisor.hh
// John McGhee
// Thu Apr  2 14:06:18 1998
//---------------------------------------------------------------------------//
// @> Defines the base class time-step advisor.
//---------------------------------------------------------------------------//

#ifndef __timestep_ts_advisor_hh__
#define __timestep_ts_advisor_hh__


#include <limits>

#include <string>

namespace rtt_timestep {

// FORWARD REFERENCES

class ts_manager; 

//===========================================================================//
// class ts_advisor - This is the base class time-step advisor.
// 
//===========================================================================//

class ts_advisor {

// NESTED CLASSES AND TYPEDEFS

  public:

// Flag to determine how the recommended timestep is to be
// used. The recommended value "dt_rec" is to be considered:

    enum usage_flag {

	min , // a lower limit
	max , // a upper limit
	req , // a required value
        last_usage  // dummy to mark end of list
    };

// DATA

  private:

    std::string name;                 //ID string
    usage_flag usage;                 //how to use dt_rec 
    bool   active;                    //on-off switch

// STATIC CLASS METHODS

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

    static std::string usage_flag_name(const int i) 
    {
	static const std::string usage_flag_names [last_usage] =
	{	"minimum",
		"maximum",
		"required"};
	return usage_flag_names[i];
    };


// CREATORS

    ts_advisor(const std::string &name_  = std::string("Unlabeled"),
	       const usage_flag usage_ = max,
	       const bool active_ = true);
    
    virtual ~ts_advisor();



//MANIPULATORS


// Turn the advisor on and off

    void activate()
    {
	active = true;
    }

    void deactivate()
    {
	active = false;
    }

// ACCESSORS

// Determine if the advisor is active or not

    bool is_active() const
    {
	return active;
    }

// Update and/or produce the recommended time-step

    virtual double get_dt_rec(const ts_manager &tsm) const = 0;

// Determine if the advisor is fit to use in
// a time-step calculation

    virtual bool advisor_usable(const ts_manager &tsm) const
    {
	return (active == true);
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
    
// Vomit the entire state of the advisor

    virtual void print_state() const = 0;

// Invariant function

    virtual bool invariant_satisfied() const
    {
	return name.length() != 0 &&
	    0      <= usage &&
	    usage  <  last_usage;
    }

// Print out advisor data

    virtual void print(const ts_manager &tsm, 
		       const bool controlling = false) const;

};

} // end of rtt_timestep namespace

#endif                          // __timestep_ts_advisor_hh__

//---------------------------------------------------------------------------//
//                              end of timestep/ts_advisor.hh
//---------------------------------------------------------------------------//
