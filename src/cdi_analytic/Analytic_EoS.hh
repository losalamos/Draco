//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/Analytic_EoS.hh
 * \author Thomas M. Evans
 * \date   Tue Oct  2 16:22:32 2001
 * \brief  Analytic_EoS class definition.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_analytic_Analytic_EoS_hh__
#define __cdi_analytic_Analytic_EoS_hh__

#include "Analytic_Models.hh"
#include "cdi/EoS.hh"
#include "ds++/SP.hh"
#include <vector>

namespace rtt_cdi_analytic
{
 
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
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Analytic_EoS : public rtt_cdi::EoS
{
  public:
    // Useful typedefs.
    typedef rtt_dsxx::SP<Analytic_EoS_Model> SP_Analytic_Model;
    typedef std::vector<double>              sf_double;

  private:
    // Analytic EoS model.
    SP_Analytic_Model analytic_model;
    
  public:
    // Constructor.
    Analytic_EoS(SP_Analytic_Model);

    // >>> INTERFACE SPECIFIED BY rtt_cdi::EoS

    // Get electron internal energy.
    double    getSpecificElectronInternalEnergy(double, double) const; 
    sf_double getSpecificElectronInternalEnergy(const sf_double &,
						const sf_double &) const;

    // Get ion internal energy.
    double    getSpecificIonInternalEnergy(double, double) const;
    sf_double getSpecificIonInternalEnergy(const sf_double &, 
					   const sf_double &) const;

    // Get electron heat capacity.
    double    getElectronHeatCapacity(double, double) const;
    sf_double getElectronHeatCapacity(const sf_double &, 
				      const sf_double &) const;

    // Get ion heat capacity.
    double    getIonHeatCapacity(double, double) const;
    sf_double getIonHeatCapacity(const sf_double &, 
				 const sf_double &) const;

    // Get the number of free electrons per ion.
    double    getNumFreeElectronsPerIon(double, double) const;
    sf_double getNumFreeElectronsPerIon(const sf_double &, 
					const sf_double &) const;

    // Get the electron thermal conductivity.
    double    getElectronThermalConductivity(double, double) const;
    sf_double getElectronThermalConductivity(const sf_double &, 
					     const sf_double &) const;
};

} // end namespace rtt_cdi_analytic

#endif                          // __cdi_analytic_Analytic_EoS_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_analytic/Analytic_EoS.hh
//---------------------------------------------------------------------------//
