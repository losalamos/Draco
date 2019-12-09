//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_EICoupling.hh
 * \author Mathew Cleveland
 * \date   March 2019
 * \brief  Analytic_EICoupling class definition.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_Analytic_EICoupling_hh__
#define __cdi_analytic_Analytic_EICoupling_hh__

#include "Analytic_Models.hh"
#include "cdi/EICoupling.hh"
#include <memory>

namespace rtt_cdi_analytic {

//===========================================================================//
/*!
 * \class Analytic_EICoupling
 *
 * \brief Derived rtt_cdi::EICoupling class for analytic electron-ion coupling
 * data.
 *
 * The Analytic_EICoupling class is a derived rtt_cdi::EICoupling class. It
 * provides analytic electron-ion coupling data.  The specific analytic
 * EICoupling model is derived from the
 * rtt_cdi_analytic::Analytic_EICoupling_Model base class.  Several pre-built
 * derived classes are provided in Analytic_Models.hh.
 *
 * Clients of this class can provide any analytic model class as long as it
 * conforms to the rtt_cdi_analytic::Analytic_EICoupling_Model interface.
 *
 * See the member functions for details about the data types and units.
 *
 * \example cdi_analytic/test/tstAnalytic_EICoupling.cc
 *
 * Example usage of Analytic_EICoupling, Analytic_EICoupling_Model, and their
 * incorporation into rtt_cdi::CDI.
 */
//===========================================================================//

class Analytic_EICoupling : public rtt_cdi::EICoupling {
public:
  // Useful typedefs.
  typedef std::shared_ptr<Analytic_EICoupling_Model> SP_Analytic_Model;
  typedef std::shared_ptr<const Analytic_EICoupling_Model> const_SP_Model;
  typedef std::vector<double> sf_double;
  typedef std::vector<char> sf_char;

private:
  // Analytic EICoupling model.
  SP_Analytic_Model analytic_model;

public:
  // Constructor.
  explicit Analytic_EICoupling(SP_Analytic_Model model_in);

  // Unpacking constructor.
  explicit Analytic_EICoupling(const sf_char &);

  // >>> ACCESSORS
  const_SP_Model get_Analytic_Model() const { return analytic_model; }

  // >>> INTERFACE SPECIFIED BY rtt_cdi::EICoupling

  // Get electron ion coupling.
  double getElectronIonCoupling(const double eTemperature,
                                const double iTemperature, const double density,
                                const double w_e, const double w_i) const;

  sf_double getElectronIonCoupling(const sf_double &vetemperature,
                                   const sf_double &vitemperature,
                                   const sf_double &vdensity,
                                   const sf_double &vw_e,
                                   const sf_double &vw_i) const;

  // Pack the Analytic_EICoupling into a character string.
  sf_char pack() const;
};

} // end namespace rtt_cdi_analytic

#endif // __cdi_analytic_Analytic_EICoupling_hh__

//---------------------------------------------------------------------------//
// end of cdi_analytic/Analytic_EICoupling.hh
//---------------------------------------------------------------------------//
