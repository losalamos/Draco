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
