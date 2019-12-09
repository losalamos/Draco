//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_CP_Eloss.cc
 * \author Kendra P. Long
 * \date   Fri Aug  2 14:28:08 2019
 * \brief  Analytic_CP_Eloss member definitions.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Analytic_CP_Eloss.hh"

namespace rtt_cdi_analytic {

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for an analytic gray eloss model.
 *
 * This constructor builds an eloss model defined by the
 * rtt_cdi_analytic::Analytic_Eloss_Model derived class argument.
 *
 * The reaction type for this instance of the class is determined by the
 * rtt_cdi::Reaction argument.
 *
 * \param model_in shared_ptr to a derived
 *                 rtt_cdi_analytic::Analytic_Eloss_Model object
 * \param target_zaid_in int32_t target particle zaid
 * \param projectile_zaid_in int32_t transporting particle zaid
 */
Analytic_CP_Eloss::Analytic_CP_Eloss(SP_Analytic_Model model_in,
                                     int32_t target_zaid_in,
                                     int32_t projectile_zaid_in)
    : analytic_model(model_in), projectile_zaid(projectile_zaid_in),
      target_zaid(target_zaid_in) {
  Ensure(analytic_model);
}

//---------------------------------------------------------------------------//
// OPACITY INTERFACE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Return a scalar eloss given a scalar temperature, density, and 
 *        particle speed.
 *
 * Given a scalar temperature/density/speed, return an eloss for the reaction
 * type specified by the constructor.  The analytic eloss model is
 * specified in the constructor (Analytic_CP_Eloss()).
 *
 * \param temperature material temperature in keV
 * \param density material density in g/cm^3
 * \param v0 incident particle speed in cm/shk
 * \return eloss (time coefficient) in shk^-1
 *
 */
double Analytic_CP_Eloss::getEloss(double temperature, double density,
                                   double v0) const {
  Require(temperature >= 0.0);
  Require(density >= 0.0);
  Require(v0 >= 0.0);

  // define return eloss
  double eloss = analytic_model->calculate_eloss(temperature, density, v0);

  Ensure(eloss >= 0.0);
  return eloss;
}

} // namespace rtt_cdi_analytic
