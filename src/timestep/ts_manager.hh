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
// recommended timesteps of its component adsisors (i.e. electron energy,
// radiation energy, ion energy, max allowed change, etc...). 
//
// Intended usage is as follows:
//
// 1) Host-code instantiates the desired time-step advisor(s).
// 2) Host-code instantiates a manager, and loads references to
///   advisor(s) into the manager list.
// ***Begin loop over time ***
//    3) Host-code responds to any user query about advisor(s), and
//       sets advisor control parameters, based on user interactive
//       input or some pre-determined script.
//    4) Host-code passes in advisor(s) ownned by each package.
//       Each package is responsible for updating recommended time step
//       of the advisors owned by that package.
//    5) Host-code updates recommended time step(s) of advisor(s) 
//       not owned by any particular package.
//    6) Manager evaluates recommendations from each advisor on the
//       list, calculates new time-step, and prints results
//       if desired.
// *** End loop over time ***
// 7) Host-code destructs manager
// 8) Host-code destructs advisor(s)
//===========================================================================//

#include <list>

#include "timestep/ts_advisor.hh"

#include "ds++/SP.hh"

#include <string>

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
    double compute_new_timestep(double dt_, int cycle_,
				double time_);

// ACCESSORS

    void print_advisors() const;
    void print_summary() const;

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

  private:

    bool invariant_satisfied() const;


};

#endif                          // __timestep_ts_manager_hh__

//---------------------------------------------------------------------------//
//                              end of timestep/ts_manager.hh
//---------------------------------------------------------------------------//
