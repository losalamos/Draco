//----------------------------------*-C++-*----------------------------------//
/*! \file   PhysicalConstants.cc
 *  \author Randy M. Roberts
 *  \brief  Provide a single place where physical constants (pi, speed of
 *          light, etc) are defined.
 *  \date   Thu Mar 19 10:04:52 1998
 *  \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "PhysicalConstants.hh"

namespace rtt_units
{

//! \brief default constructor use SI units (kg, m, seconds, degree K)
PhysicalConstants::PhysicalConstants()
    : d_planck         ( planckSI ),
      d_gasConstant    ( gasConstantSI ),
      d_boltzmann      ( boltzmannSI ),
      d_electronCharge ( electronChargeSI ),
      d_cLight         ( cLightSI ),
      d_stefanBoltzmann( stefanBoltzmannSI ),
      d_gravitationalConstant( gravitationalConstantSI ),
      d_accelerationFromGravity( accelerationFromGravitySI ),
      d_faradayConstant( faradayConstantSI ),
      d_permeabilityOfVacuum( permeabilityOfVacuumSI ),
      d_permittivityOfFreeSpace( permittivityOfFreeSpaceSI ),
      d_electronMass   ( electronMassSI ),
      d_protonMass     ( protonMassSI )
{
    // empty
}

//! \brief constuctor based on rtt_units::UnitSystem
PhysicalConstants::PhysicalConstants( UnitSystem const & u )
    : d_planck          ( planckSI  * u.e() * u.t() ),
      d_gasConstant     ( gasConstantSI * u.e() / u.T() ),
      d_boltzmann       ( d_gasConstant / AVOGADRO ),
      d_electronCharge  ( electronChargeSI ),
      d_cLight          ( cLightSI * u.v() ),
      d_stefanBoltzmann ( 
	  stefanBoltzmannSI * u.p() / pow(u.L(),2) / pow(u.T(),4) ),
      d_gravitationalConstant( 
	  gravitationalConstantSI * u.f() * pow(u.L(),2) / pow(u.M(),2) ),
      d_accelerationFromGravity( accelerationFromGravitySI * u.a() ),
      d_faradayConstant( AVOGADRO * d_electronCharge ),
      d_permeabilityOfVacuum( permeabilityOfVacuumSI / u.L() ),
      d_permittivityOfFreeSpace( 
	  1.0 / d_permeabilityOfVacuum / pow(d_cLight,2) ),
      d_electronMass    ( electronMassSI * u.M() ),
      d_protonMass      ( protonMassSI   * u.M() )
{
    // empty
}

} // end namespace rtt_units

//---------------------------------------------------------------------------//
//                  end of units/PhysicalConstants.cc
//---------------------------------------------------------------------------//
