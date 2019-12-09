//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   dummy_package.hh
 * \author John McGhee
 * \date   Thu Aug 27 07:48:41 1998
 * \brief  A dummy package to exercize the field time-step advisors.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __timestep_test_dummy_package_hh__
#define __timestep_test_dummy_package_hh__

#include "ds++/config.h"
#include <memory>
#include <vector>

// FORWARD REFERENCES

namespace rtt_timestep {
class ts_manager;
class field_ts_advisor;
} // namespace rtt_timestep

namespace rtt_timestep_test {

//===========================================================================//
/*! 
 * \class dummy_package 
 * \brief Exercizes the field time-step advisors.
 *
 * This class serves as an example of how any particular package can make use
 * of the time-step manaager/advisor utility. It also exercizes the field
 * time-step advisors.
 */
//===========================================================================//

class DLL_PUBLIC_timestep_test dummy_package {
public:
  // CREATORS

  explicit dummy_package(rtt_timestep::ts_manager &tsm_);
  ~dummy_package(void);

  // MANIPULATORS

  //! \brief Method to advance the time-step
  void advance_state(void);

  //! \brief Provide recommended dt values (TE).
  double get_dt_rec_te(void) const;

  //! \brief Provide recommended dt values (TI).
  double get_dt_rec_ti(void) const;

  //! \brief Provide recommended dt values (RI).
  double get_dt_rec_ri(void) const;

private:
  std::vector<double> element_wise_multiply(double const a,
                                            std::vector<double> const &v);

  // DATA

  rtt_timestep::ts_manager &tsm;
  std::shared_ptr<rtt_timestep::field_ts_advisor> sp_te;
  std::shared_ptr<rtt_timestep::field_ts_advisor> sp_ti;
  std::shared_ptr<rtt_timestep::field_ts_advisor> sp_ri;
};

} // end namespace rtt_timestep_test

#endif // __timestep_test_dummy_package_hh__

//---------------------------------------------------------------------------//
// end of dummy_package.hh
//---------------------------------------------------------------------------//
