//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   timestep/field_ts_advisor.hh
 * \author John McGhee
 * \date   Thu Apr  2 14:06:18 1998
 * \brief  Header file for the field time-step advisor class.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __timestep_field_ts_advisor_hh__
#define __timestep_field_ts_advisor_hh__

#include "ts_advisor.hh"

namespace rtt_timestep {

//===========================================================================//
/*!
 * \brief Estimates a new timestep based on current fields.
 *
 * \sa The ts_manager class provides a means to manage this advisor.
 * The \ref overview_timestep page gives a summary of the Draco time 
 * step control utilities. 
 *
 * This class provides a means to estimate a
 * new timestep based on the current dQ/dt where Q is some field quantity
 * to be controlled (i.e. temperature, energy, particle number density,
 * etc....), and dt is the current timestep. 
 *
 * The recommended timestep will be
 * based on a norm of a control function "alpha", where
 * alpha = abs(del_Q/Q_old), and where Q is a field  being monitored,
 * i.e. temperature, energy, particles, etc. 
 * Alpha is computed point-wise in the field, then a vector
 * norm is applied to the point-wise values. The type of norm used is selected
 * by the update method flag.
 */
//===========================================================================//
class DLL_PUBLIC_timestep field_ts_advisor : public ts_advisor {

  // NESTED CLASSES AND TYPEDEFS

public:
  //! Selects the vector norm used to evaluate changes in the control field.
  enum update_method_flag {
    inf_norm, //!< infinity norm (max)
    a_mean,   //!< arithmetic mean
    q_mean,   //!< Q weighted mean
    rc_mean,  //!< relative change (alpha) weighted mean
    rcq_mean, //!< product of Q and alpha weighted mean
    last_umf  //!< dummy to mark end of list
  };

  // DATA

private:
  update_method_flag update_method; //!< update method for dt_rec
  double fc_value;                  //!< frac change  value for field advisor
  double floor_value;               //!< floor value for field advisor
  int cycle_at_last_update;         //!< problem time-cycle index at last update
  double dt_rec;                    //!< the recommended time-step

  // STATIC CLASS METHODS

public:
  //! Returns the name of the update method flag requested
  /*! \param i the update_method_flag for which a name is desired
   *  \return the name of the vector norm associated with i 
   */
  static std::string update_method_flag_name(const int i) {
    static const std::string update_method_flag_names[last_umf] = {
        "infinity norm", "arithmetic mean", "weighted by field value",
        "weighted by relative change", "field value and relative change"};
    return update_method_flag_names[i];
  };

  // CREATORS

  //! Constructs a field time step advisor
  /*! \param name_ A unique name for the advisor
   *  \param usage_  Specifies how the advisor is to be used
   *  \param update_method_ Specifies the update method to be used
   *  \param fc_value_ Specifies the desired fractional change
   *  \param floor_value_ Specifies the floor value
   *  \param active_ Turns the advisor on/off
   */
  field_ts_advisor(const std::string &name_ = std::string("Unlabeled"),
                   const usage_flag usage_ = max,
                   const update_method_flag update_method_ = inf_norm,
                   const double fc_value_ = 0.1,
                   const double floor_value_ = ts_small(),
                   const bool active_ = true);

  //! Destroys a field time step advisor
  ~field_ts_advisor(void){/*empty*/};

  // MANIPULATORS

  /*! \brief The floor is computed as the max of frac*y1
   *
   * A utility function to set the floor value as a fraction of the max
   * value in a field 
   *
   *  \param y1 The field to be examined
   *  \param frac The fractional value to be applied to the field
   */
  template <class FT> void set_floor(FT const &y1, double frac = 0.001);

  //! Calculates a new recommended time-step for a field advisor.
  /*! \param tsm the time step manager in which the advisor resides 
      \param q_old  the field value at the beginning of the current time-step,
      \param q_new the field value at the end of the current time-step.
  */
  template <class FT>
  void update_tstep(const ts_manager &tsm, const FT &q_old, const FT &q_new);

  //! Sets the fractional change value
  /*! \param value_ the value to set the fractional change
     */
  void set_fc(const double value_ = 0.1) { fc_value = value_; }

  //! Sets the floor value
  /*! \param value_ the value to set the floor
   */
  void set_floor(const double value_ = ts_small()) { floor_value = value_; }

  //! Sets the update method
  /*! \param flag  the flag to set
   */
  void set_update_method(const update_method_flag flag = inf_norm) {
    update_method = flag;
  }

  // ACCESSORS

  //! Produces the recommended time-step
  /*! \param tsm the time step manager in which the advisor resides 
   *  \return the recommended timestep
   */
  double get_dt_rec(const ts_manager &tsm) const;

  //! Determines if the advisor is fit to use in a time-step calculation
  /*! \param tsm the time step manager in which the advisor resides 
   *  \return true if the advisor is usable
   */
  bool advisor_usable(const ts_manager &tsm) const;

  //! Prints advisor state
  /*! \return prints the advisor internal state to std out
   */
  void print_state(std::ostream &out = std::cout) const;

  //! Invariant function
  /*! \return true if the invariant is satisfied.
   */
  bool invariant_satisfied() const;
};

} //end namespace rtt_timestep

#endif // __timestep_field_ts_advisor_hh__

//---------------------------------------------------------------------------//
// end of field_ts_advisor.hh
//---------------------------------------------------------------------------//
