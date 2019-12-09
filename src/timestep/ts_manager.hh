//----------------------------------*-C++-*----------------------------------//
/*!
 * \file    timestep/ts_manager.hh
 * \author  John McGhee
 * \date    Mon Apr  6 17:22:53 1998
 * \brief   Header file for the manager utility for time-step advisors.
 * \note    Copyright (C) 2016-2019 Triad National Security, LLC.
 *          All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __timestep_ts_manager_hh__
#define __timestep_ts_manager_hh__

#include "ts_advisor.hh"
#include <list>
#include <memory>

namespace rtt_timestep {

//===========================================================================//
/*!
 * \brief Manages a list of time-step advisors.
 *
 * \sa The ts_advisor class provides the advisors to be registerd the the
 *      ts_manager class. Also, the \ref overview_timestep page provides useful
 *      info.
 *
 * Calculates a new timestep based on the recommended time-steps of its
 * component advisors (i.e. electron energy, radiation energy, ion energy, max
 * allowed change, etc...).
 */
//===========================================================================//
class DLL_PUBLIC_timestep ts_manager {

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
  int cycle;
  //! name of the advisor in control
  std::string controlling_advisor;
  //! a list of Smart Pointers to time step advisors
  std::list<std::shared_ptr<ts_advisor>> advisors;

  // CREATORS

public:
  //! Creates a timestep manager
  ts_manager();

  //! Destroys a timestep manager
  ~ts_manager(){/* empty */};

  // MANIPULATORS

  /*!
   * \brief Adds a timestep advisor to a RTT timestep manager
   * \param new_advisor the new advisor to be added
   */
  void add_advisor(const std::shared_ptr<ts_advisor> &new_advisor);

  /*!
   * \brief Removes a timestep advisor from a RTT timestep manager
   * \param advisor_to_remove the advisor to be removed
   */
  void remove_advisor(const std::shared_ptr<ts_advisor> &advisor_to_remove);

  /*!
   * \brief Sets timestep, cycle number, and problem time
   * \param dt_ the timestep (time)
   * \param cycle_ the time cycle
   * \param time_ the problem time (time)
   */
  void set_cycle_data(double dt_, int cycle_, double time_);

  /*!
   * \brief Computes a new timestep based on the recommendations of each advisor
   * \return the recommended timestep
   */
  double compute_new_timestep();

  // ACCESSORS

  /*!
   * \brief Prints advisor names
   * \return prints a list of the advisor names to std out
   */
  void print_advisors() const;

  /*!
   * \brief Prints a concise summary of the manager status
   * \return prints a summary to std out
   */
  void print_summary() const;

  /*!
   * \brief Prints advisor information
   * Prints a detailed listing of each advisors internal state to std out
   */
  void print_adv_states() const;

  /*!
   * \brief Returns the recommendation for the next timestep
   * \return the recommended timestep
   */
  double get_dt_new() const { return dt_new; }

  //! Returns the current problem time
  double get_time() const { return time; }

  //! Returns the current timestep
  double get_dt() const { return dt; }

  //! Returns the current time cycle number
  int get_cycle() const { return cycle; }

  /*!
   * \brief Returns the controlling advisor
   * \return The name of the controlling advisor
   */
  std::string get_controlling_advisor() const { return controlling_advisor; }

  /*!
   * \brief Defines the timestep manager invariant function
   * \return True if the invariant is satisfied
   */
  bool invariant_satisfied() const;
};

} // namespace rtt_timestep

#endif // __timestep_ts_manager_hh__

//---------------------------------------------------------------------------//
// end of ts_manager.hh
//---------------------------------------------------------------------------//
