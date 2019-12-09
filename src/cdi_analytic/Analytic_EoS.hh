//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_EoS.hh
 * \author Thomas M. Evans
 * \date   Tue Oct  2 16:22:32 2001
 * \brief  Analytic_EoS class definition.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_Analytic_EoS_hh__
#define __cdi_analytic_Analytic_EoS_hh__

#include "Analytic_Models.hh"
#include "cdi/EoS.hh"
#include <memory>

namespace rtt_cdi_analytic {

//===========================================================================//
/*!
 * \class Analytic_EoS
 *
 * \brief Derived rtt_cdi::EoS class for analytic Equation of State data.
 *
 * The Analytic_EoS class is a derived rtt_cdi::EoS class. It provides
 * analytic Equation of State (EoS) data.  The specific analytic EoS model is
 * derived from the rtt_cdi_analytic::Analytic_EoS_Model base class.  Several
 * pre-built derived classes are provided in Analytic_Models.hh.
 *
 * Clients of this class can provide any analytic model class as long as it
 * conforms to the rtt_cdi_analytic::Analytic_EoS_Model interface.
 *
 * See the member functions for details about the data types and units.
 *
 * \example cdi_analytic/test/tstAnalytic_EoS.cc
 *
 * Example usage of Analytic_EoS, Analytic_EoS_Model, and their incorporation
 * into rtt_cdi::CDI.
 */
// revision history:
// -----------------
// 0) original
//
//===========================================================================//

class DLL_PUBLIC_cdi_analytic Analytic_EoS : public rtt_cdi::EoS {
public:
  // Useful typedefs.
  typedef std::shared_ptr<Analytic_EoS_Model> SP_Analytic_Model;
  typedef std::shared_ptr<const Analytic_EoS_Model> const_SP_Model;
  typedef std::vector<double> sf_double;
  typedef std::vector<char> sf_char;

private:
  // Analytic EoS model.
  SP_Analytic_Model analytic_model;

public:
  // Constructor.
  explicit Analytic_EoS(SP_Analytic_Model model_in);

  // Unpacking constructor.
  explicit Analytic_EoS(const sf_char &);

  // >>> ACCESSORS
  const_SP_Model get_Analytic_Model() const { return analytic_model; }

  // >>> INTERFACE SPECIFIED BY rtt_cdi::EoS

  // Get electron internal energy.
  double getSpecificElectronInternalEnergy(double, double) const;
  sf_double getSpecificElectronInternalEnergy(const sf_double &,
                                              const sf_double &) const;

  // Get ion internal energy.
  double getSpecificIonInternalEnergy(double, double) const;
  sf_double getSpecificIonInternalEnergy(const sf_double &,
                                         const sf_double &) const;

  // Get electron heat capacity.
  double getElectronHeatCapacity(double, double) const;
  sf_double getElectronHeatCapacity(const sf_double &, const sf_double &) const;

  // Get ion heat capacity.
  double getIonHeatCapacity(double, double) const;
  sf_double getIonHeatCapacity(const sf_double &, const sf_double &) const;

  // Get the number of free electrons per ion.
  double getNumFreeElectronsPerIon(double, double) const;
  sf_double getNumFreeElectronsPerIon(const sf_double &,
                                      const sf_double &) const;

  // Get the electron thermal conductivity.
  double getElectronThermalConductivity(double, double) const;
  sf_double getElectronThermalConductivity(const sf_double &,
                                           const sf_double &) const;

  // Get the new Te, given delta Ue, Te0.
  double getElectronTemperature(double /*rho*/, double Ue,
                                double Tguess = 1.0) const;

  // Get the new Ti, given delta Uic, Ti0.
  double getIonTemperature(double /*rho*/, double Uic,
                           double Tguess = 1.0) const;

  // Pack the Analytic_EoS into a character string.
  sf_char pack() const;
};

} // end namespace rtt_cdi_analytic

#endif // __cdi_analytic_Analytic_EoS_hh__

//---------------------------------------------------------------------------//
// end of cdi_analytic/Analytic_EoS.hh
//---------------------------------------------------------------------------//
