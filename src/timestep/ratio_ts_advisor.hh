//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   timestep/ratio_ts_advisor.hh
 * \brief  Header file for the ratio time-step advisor class.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __timestep_ratio_ts_advisor_hh__
#define __timestep_ratio_ts_advisor_hh__

#include "ts_advisor.hh"
#include "ds++/Assert.hh"

namespace rtt_timestep {

//===========================================================================//
/*! 
 * \brief Calculates a new timestep as a ratio of the current time-step
 *
 * \sa The ts_manager class provides a means to manage this advisor.
 * The \ref overview_timestep page gives a summary of the Draco time
 * step control utilities. 
 *
 * This class provides a means to calculate a
 * new timestep as a ratio of the current time-step. This is useful to
 * limit the rate of change in the time-step from one time cycle to the 
 * next. The recommendation for the new time step is computed as 
 * current_dt*ratio. 
 */
//===========================================================================//
class DLL_PUBLIC_timestep ratio_ts_advisor : public ts_advisor {

  // DATA

private:
  //! Value used to compute ratio advisor
  double ratio_value;

  // CREATORS
public:
  //! Constructs a ratio time step advisor
  /*! \param name_ a unique name for the advisor
   *  \param usage_ how the advisor is to be used
   *  \param ratio_value_ the value of the ratio to be used
   *  \param active_ turns the advisor on/off
   */
  ratio_ts_advisor(const std::string &name_ = std::string("Unlabeled"),
                   const usage_flag usage_ = max,
                   const double ratio_value_ = 1.20, const bool active_ = true);

  //! Destroys a ratio time step advisor
  ~ratio_ts_advisor();

  // MANIPULATORS

  //! Set the ratio value
  /*! \param value_ the value of the desired ratio, (value > 0.)
   */
  void set_ratio(const double value_ = 1.2) {
    Require(value_ > 0.0);
    ratio_value = value_;
  }

  // ACCESSORS

  //! Returns the recommended time-step
  /*! \return the time step recommended by this advisor  
   *  \param tsm the time step manager in which the advisor resides
   */
  double get_dt_rec(const ts_manager &tsm) const;

  //! Prints state
  /*! \return Prints the internal state of the advisor to std out 
   */
  void print_state(std::ostream &out = std::cout) const;

  //! Invariant function
  /*! \return True if the invariant is satisfied 
   */
  bool invariant_satisfied() const;

  //! Returns the current ratio_value;
  double get_ratio() const {
    Ensure(ratio_value > 0.0);
    return ratio_value;
  }
};

} // namespace rtt_timestep

#endif // __timestep_ratio_ts_advisor_hh__

//---------------------------------------------------------------------------//
// end of ratio_ts_advisor.hh
//---------------------------------------------------------------------------//
