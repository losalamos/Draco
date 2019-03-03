//----------------------------------*-C++-*----------------------------------//
/*!
 * \file timestep/target_ts_advisor.hh
 * \brief Header file for the target time-step advisor class.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __timestep_target_ts_advisor_hh__
#define __timestep_target_ts_advisor_hh__

#include "ts_advisor.hh"

namespace rtt_timestep {

//===========================================================================//
/*!
 * \brief  Calculates a new time-step required to hit some "target" 
 * problem time.
 *
 * \sa The ts_manager class provides a means to manage this advisor.
 * The \ref overview_timestep page gives a summary of the Draco time 
 * step control utilities. 
 *
 * This class provides a means to calculate 
 * a new time-step required to hit some "target" problem time. This
 * is useful to assure that graphics dumps/IO/etc. occur precisely at
 * a predetermined problem time. 
 */
//===========================================================================//
class DLL_PUBLIC_timestep target_ts_advisor : public ts_advisor {

  // DATA

private:
  //! Target time (time)
  double target_value;

  // CREATORS

public:
  //! Constructs a target time step advisor
  /*! \param name_ a unique name for the advisor
   *  \param usage_ how the advisor is to be used
   *  \param target_value_ the problem target time (time)
   *  \param active_ turns the advisor on/off
   */
  target_ts_advisor(const std::string &name_ = std::string("Unlabeled"),
                    const usage_flag usage_ = max,
                    const double target_value_ = -large(),
                    const bool active_ = true);

  //! Destroys a target time step advisor
  ~target_ts_advisor(void){/*empty*/};

  // MANIPULATORS

  //! Set the target value
  /*! \param value_ the target value (time)
   */
  void set_target(double const value_ = -large()) { target_value = value_; }

  // ACCESSORS

  //! Get the target value
  double get_target() { return target_value; }

  //! Returns the recommended time-step
  /*! \param tsm the time step manager in which the advisor resides 
   *  \return the time step recommended by this advisor 
   */
  double get_dt_rec(const ts_manager &tsm) const;

  //! Prints state
  /*! Prints the internal state of the advisor to std out
   */
  void print_state(std::ostream &out = std::cout) const;

  //! Invariant function
  /*! \return True if the invariant is satisfied
   */
  bool invariant_satisfied() const;
};

} // namespace rtt_timestep

#endif // __target_timestep_ts_advisor_hh__

//---------------------------------------------------------------------------//
// end of target_ts_advisor.hh
//---------------------------------------------------------------------------//
