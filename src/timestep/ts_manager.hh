//----------------------------------*-C++-*----------------------------------//
// ts_manager.hh
// John McGhee
// Mon Apr  6 17:22:53 1998
//---------------------------------------------------------------------------//
// @> Defines a manager utility for time-step advisors.
//---------------------------------------------------------------------------//

/*! 
 * \file
 * \brief Defines a manager utility for time-step advisors.
 *
 * \author <a href="http://www.lanl.gov/home/mcghee">
 *  John McGhee</a>
 *
 * \date Mon Apr  6 17:22:53 1998
 *
 */

/*!
 * \example main.cc
 * The following code provides an example of how to use the timestep manager
 * utility.
 * \include test_timestep.hh
 * \include test_timestep.cc
 */

/*!
 * \example dummy_package.cc
 * The following code provides a dummy package for use with the
 * test_timestep example.
 * \include dummy_package.hh
 * \include test_timestep_pt.cc
 */

/*!
 * \example test_timestep.out
 * The following code is a sample output from the test_timestep example.
 * It contains representative output from most of the printing and summary
 * I/O utilities.
 */

/*!
 * \page overview Overview of the Draco Time Step Manager
 *
 * <h3> Introduction </h3>
 * The classes contained in the rtt_timestep name space are designed
 * to provide automatic control of the time step size in a time dependent
 * numerical simulation. In a typical simulation, several physics
 * packages may be repeatedly called from a host code to advance the
 * overall system state. 
 * A time step manager class (rtt_timestep::ts_manager) is provided 
 * to contain any number of
 * time step advisors (rtt_timestep::ts_advisor) from any package 
 * in the simulation. Each advisor
 * provides a recommended time step based on some physical parameter
 * specific to the package. After collecting information from each
 * advisor, the time step manager considers all the inputs and provides
 * an overall recommendation for the size of the time step for the
 * next time cycle.  Capabilites are provided to:
 * <ul>
 *  <li> limit the min and max timestep,
 *  <li> force the time step to a user input value,
 *  <li> display informational advisors,
 *  <li> print a variety of summaries,
 *  <li> limit the rate of change in the time step
 *  <li> control the time step based on estimated rate of change in
 *       a user selected physical parameter,
 *  <li> attempt to hit a "target" time. 
 * </ul> 
 *
 * <h3> Intended Usage </h3>
 * The intended usage is as follows:
 *
 * 1) Host-code instantiates the desired global (i.e. package
 *    independent) time-step advisor(s).
 *
 * 2) Host-code instantiates a ts_manager, and loads references to
 *    advisor(s) created in 1) above into the manager list.
 *
 * 3) Each package in the host code, at package construction, 
 *    is passed in the manager, and instantiates any 
 *    desired time-step advisors and loads references 
 *    to the new advisors into the time-step manager.
 *
 * ***Begin loop over time ***
 *
 *    4) Host-code responds to any user query about advisor(s), and
 *       sets advisor control parameters, based on user interactive
 *       input or some pre-determined script.
 *
 *    5) Host-code sets cycle data of manager (dt, cycle#, time).
 *
 *    6) Each package  updates, if required, the recommended time
 *       step for each advisor "owned" by that package.
 *
 *    7) Host-code updates recommended time step(s) of global 
 *       advisor(s) not owned by any particular package, if required.
 *
 *    8) Manager evaluates recommendations from each advisor on the
 *       list, calculates new time-step, and prints results
 *       if desired.
 *
 * *** End loop over time ***
 *
 * 9)  At package destruction, each package removes the advisors
 *     created in 3) above from the manager.
 *
 * 10) Host-code removes global advisors from the manager list.
 *
 * 11) Host-code destructs manager.
 *
 * <h3> Other Draco Packages </h3>
 * The time step manager utility uses the Draco ds++ services library, 
 * and the Draco C4 communication library.
 */

#ifndef __timestep_ts_manager_hh__
#define __timestep_ts_manager_hh__

#include <list>

#include "ts_advisor.hh"

#include "ds++/SP.hh"

#include <string>

//! RTT time step namespace
/*!
 * Provides namspace protection for the Draco (RTT) time step 
 * control utilities.
 *\sa The ts_manager and ts_advisor classes provide most of the
 *    functionality of the namspace. The \ref overview page presents
 *    a summary of the capabilities provided within the namespace.
 */
namespace rtt_timestep {

//===========================================================================//
//! Manages a list of time-step advisors.
/*!
 * \sa  The ts_advisor class provides the advisors to be registerd
 *      the the ts_manager class. Also, the \ref overview page provides
 *      useful info.
 *
 * Calculates a new timestep based on the 
 * recommended time-steps of its component advisors (i.e. electron energy,
 * radiation energy, ion energy, max allowed change, etc...). 
 */
//===========================================================================//
class ts_manager {

// NESTED CLASSES AND TYPEDEFS

// DATA

  private:

    //! the recommendation for the next time-step (time)
    double dt_new; 
    //! problem time at the end of current cycle  (time)
    double time;  
    //! the current time-step (time) 
    double dt;  
    //! current cycle number   
    int    cycle;  
    //! name of the advisor in control
    std::string controlling_advisor; 
    //! a list of Smart Pointers to time step advisors
    std::list < dsxx::SP<ts_advisor> > advisors; 

// CREATORS

  public:
    //! Creates a timestep manager
    ts_manager();

    //! Destroys a timestep manager
    ~ts_manager();

// MANIPULATORS

    //! Adds a timestep advisor to a RTT timestep manager
    /*! \param new_advisor the new advisor to be added
     */
    void add_advisor(const dsxx::SP<ts_advisor> &new_advisor);

    //! Removes a timestep advisor from a RTT timestep manager
    /*! \param advisor_to_remove the advisor to be removed 
     */
    void remove_advisor(const dsxx::SP<ts_advisor> &advisor_to_remove);

    //! Sets timestep, cycle number, and problem time
    /*! \param dt_ the timestep (time)
        \param cycle_ the time cycle
        \param time_ the problem time (time)
    */
    void set_cycle_data(double dt_, int cycle_, double time_);

    //! Computes a new timestep based on the recommendations of each advisor
    /*! \return the recommended timestep
     */
    double compute_new_timestep();

// ACCESSORS

    //! Prints advisor names
    /*! \return prints a list of the advisor names to std out
     */
    void print_advisors() const;

    //! Prints a concise summary of the manager status
    /*! \return prints a summary to std out
     */
    void print_summary() const;

    //! Prints advisor information
    /*! Prints a detailed listing of each advisors internal state to std out
     */
    void print_adv_states() const;

    //! Returns the recommendation for the next timestep 
    /*! \return the recommended timestep
     */
    double get_dt_new() const
    {return dt_new;}

    //! Returns the current problem time
    double get_time() const
    {return time;}

    //! Returns the current timestep
    double get_dt() const
    {return dt;}

    //! Returns the current time cycle number
    int get_cycle() const
    {return cycle;}

    //! Returns the controlling advisor
    /*! \return The name of the controlling advisor
     */
    std::string get_controlling_advisor() const
    {return controlling_advisor;}

    //! Defines the timestep manager invariant function
    /*! \return True if the invariant is satisfied
     */
    bool invariant_satisfied() const;

};

} // end of rtt_timestep namespace

#endif                          // __timestep_ts_manager_hh__

//---------------------------------------------------------------------------//
//                              end of ts_manager.hh
//---------------------------------------------------------------------------//
