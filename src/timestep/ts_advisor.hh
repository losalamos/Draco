//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   timestep/ts_advisor.hh
 * \author John McGhee
 * \date   Thu Apr  2 14:06:18 1998
 * \brief  Header file for the base class time-step advisor.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __timestep_ts_advisor_hh__
#define __timestep_ts_advisor_hh__

#include "ds++/config.h"
#include <iostream>
#include <limits>
#include <string>

namespace rtt_timestep {

// FORWARD REFERENCES

class ts_manager;

//===========================================================================//
/*!
 * \brief Defines the base class time-step advisor.
 *
 * \sa The ts_manager class provides a means to manage this advisor. The \ref 
 *     overview_timestep page gives a summary of the Draco time step control
 *     utilities.
 */
//===========================================================================//
class DLL_PUBLIC_timestep ts_advisor {

  // NESTED CLASSES AND TYPEDEFS

public:
  /*! 
   * \brief Flag to determine how the recommended timestep is to be used.
   *
   * The recommended value "dt_rec" is to be used as indicated by the
   * enumeration value selected
   */
  enum usage_flag {
    min,       //!< use as a lower limit
    max,       //!< use as a upper limit
    req,       //!< use as a required value
    last_usage //!< dummy to mark end of list
  };

  // DATA

private:
  std::string name; //!< ID string
  usage_flag usage; //!< how to use dt_rec
  bool active;      //!< on-off switch

  // STATIC CLASS METHODS

public:
  //! Returns a number close to machine precision
  static double eps() { return 100. * std::numeric_limits<double>::epsilon(); }

  //! Returns a small number
  static double ts_small() { return 100. * std::numeric_limits<double>::min(); }

  //! Returns a large number
  static double large() { return 0.01 * std::numeric_limits<double>::max(); }

  //! Returns the name of the usage flag requested
  static std::string usage_flag_name(const int i) {
    static const std::string usage_flag_names[last_usage] = {
        "minimum", "maximum", "required"};
    return usage_flag_names[i];
  };

  // CREATORS

  /*!
   * \brief Constucts the advisor
   *
   * \param name_ the name of the advisor
   * \param usage_ Specifies how the advisor is to be used
   * \param active_ turns the advisor on/off
   */
  ts_advisor(const std::string &name_ = std::string("Unlabeled"),
             const usage_flag usage_ = max, const bool active_ = true);

  //! Destroy the advisor
  virtual ~ts_advisor(void){/*empty*/};

  //MANIPULATORS

  //! Turn the advisor on
  // void activate() { active = true; }

  //! Turn the advisor off
  void deactivate() { active = false; }

  // ACCESSORS

  //! Determine if the advisor is active or not
  bool is_active() const { return active; }

  /*!
   * \brief Update and/or produce the time-step recommended by this advisor
   * \param tsm the timestep manager in which the advisor resides
   */
  virtual double get_dt_rec(const ts_manager &tsm) const = 0;

  /*!
   * \brief Determine if the advisor is fit to use in a time-step calculation
   *
   * Derived classes will have parameter 'tsm', the timestep manager in which 
   * the advisor resides.
   */
  virtual bool advisor_usable(const ts_manager & /*tsm*/) const

  {
    return (active == true);
  }

  //! Get the usage
  usage_flag get_usage() const { return usage; }

  //! Get the name
  const std::string &get_name() const { return name; }

  //! Vomit the entire internal state of the advisor to std out
  virtual void print_state(std::ostream &out = std::cout) const = 0;

  //! True if the invariant is satisfied.
  virtual bool invariant_satisfied() const {
    return name.length() != 0 && 0 <= usage && usage < last_usage;
  }

  /*!
   * \brief Print out advisor data
   *
   * \param tsm the timestep manager in which the advisor resides
   * \param controlling flags the advisor as the controlling advisor
   */
  virtual void print(const ts_manager &tsm,
                     const bool controlling = false) const;
};

} // namespace rtt_timestep

#endif // __timestep_ts_advisor_hh__

//---------------------------------------------------------------------------//
// end of ts_advisor.hh
//---------------------------------------------------------------------------//
