//----------------------------------*-C++-*----------------------------------//
// ts_manager.hh
// John McGhee
// Mon Apr  6 17:22:53 1998
//---------------------------------------------------------------------------//
// @> Defines a manager utility for time-step advisors.
//---------------------------------------------------------------------------//

#ifndef __timestep_ts_manager_hh__
#define __timestep_ts_manager_hh__

//===========================================================================//
// class ts_manager - Manages a list of time-step advisors.
//
//
// Calculates a new timestep based on the 
// recommended time-steps of its component advisors (i.e. electron energy,
// radiation energy, ion energy, max allowed change, etc...). 
//
// Intended usage is as follows:
//
// 1) Host-code instantiates the desired global (i.e. package
//    independent) time-step advisor(s).
// 2) Host-code instantiates a manager, and loads references to
//    advisor(s) created in 1) above into the manager list.
// 3) Each package in the host code, at package construction, 
//    is passed in the manager, and instantiates any 
//    desired time-step advisors and loads references 
//    to the new advisors into the time-step manager.
// ***Begin loop over time ***
//    4) Host-code responds to any user query about advisor(s), and
//       sets advisor control parameters, based on user interactive
//       input or some pre-determined script.
//    5) Host-code sets cycle data of manager (dt, cycle#, time).
//    6) Each package  updates, if required, the recommended time
//       step for each advisor "owned" by that package.
//    7) Host-code updates recommended time step(s) of global 
//       advisor(s) not owned by any particular package, if required.
//    8) Manager evaluates recommendations from each advisor on the
//       list, calculates new time-step, and prints results
//       if desired.
// *** End loop over time ***
// 9)  At package destruction, each package removes the advisors
//     created in 3) above from the manager.
// 10) Host-code removes global advisors from the manager list.
// 11) Host-code destructs manager.
//===========================================================================//

#include <list>

#include "timestep/ts_advisor.hh"

#include "ds++/SP.hh"

#include <string>

namespace rtt_timestep {

class ts_manager {

// NESTED CLASSES AND TYPEDEFS

// DATA

  private:

    double dt_new; // the recommendation for the next time-step (time)
    double time;   // problem time at the end of current cycle  (time)
    double dt;     // the current time-step (time)
    int    cycle;  // current cycle number
    std::string controlling_advisor; // name of the advisor in control
    std::list < dsxx::SP<ts_advisor> > advisors; // a list of Smart Pointers to
// time step advisors

// CREATORS

  public:

    ts_manager();
    ~ts_manager();

// MANIPULATORS

    void add_advisor(const dsxx::SP<ts_advisor> &new_advisor);
    void remove_advisor(const dsxx::SP<ts_advisor> &advisor_to_remove);
    void set_cycle_data(double dt_, int cycle_, double time_);
    double compute_new_timestep();

// ACCESSORS

    void print_advisors() const;
    void print_summary() const;
    void print_adv_states() const;

    double get_dt_new() const
    {return dt_new;}

    double get_time() const
    {return time;}

    double get_dt() const
    {return dt;}

    int get_cycle() const
    {return cycle;}

    std::string get_controlling_advisor() const
    {return controlling_advisor;}

    bool invariant_satisfied() const;

};

} // end of rtt_timestep namespace

#endif                          // __timestep_ts_manager_hh__

//---------------------------------------------------------------------------//
//                              end of timestep/ts_manager.hh
//---------------------------------------------------------------------------//
