//----------------------------------*-C++-*----------------------------------//
/*! \file   PhysicalConstantsSI.hh
 *  \author Kelly Thompson, Kent G. Budge
 *  \brief  Provide a single place where physical constants (pi, speed of
 *          light, etc) are defined in SI units.
 *  \date   Fri Nov 07 10:04:52 2003
 *  \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *          All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: PhysicalConstants.hh 7431 2015-02-24 17:10:44Z kellyt $
//---------------------------------------------------------------------------//

#ifndef __units_PhysicalConstantsSI_hh__
#define __units_PhysicalConstantsSI_hh__

#include "MathConstants.hh"
#include "UnitSystem.hh"

//! \namespace rtt_units Namespace for units and physical constants
namespace rtt_units {
// Base physical constants in SI units:

//    m - meters, kg - kilograms, s - seconds, K - kelvin
//    W - watts, J - joules, C - coulombs, F - farads
//    mol - mole

//---------------------------------------------------------------------------//
// FUNDAMENTAL CONSTANTS
//
// NIST Reference on Constants, Units and Uncertainty
// CODATA internationally recommended values of the Fundamental Physical
// Constants, http://physics.nist.gov/cuu/Constants/
//
// The units of these fundamental constants should be factors of 10X different
// from the official NIST 2010 CODATA report data to allow for easy comparison
// between these values and the NIST data.
//
// Fundamental constants are listed first.
// Derived constants are listed second.
// Actual data is placed in a user-defined type for C-interoperatbility.
//
//---------------------------------------------------------------------------//

//! [c] SPEED OF LIGHT (M/S)
// exact value by NIST definition
static double const cLightSI = 2.99792458e8; // m s^-1

//! [Na] AVOGADRO'S NUMBER ("entities"/MOLE)
// Wikipedia (2013-12-3) == NIST Codata 2010
static double const AVOGADRO = 6.02214129e23; // entities/mol

//! [h] Planck constant ( Energy seconds )
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 4.4e-8)
static double const planckSI = 6.62606957e-34; // J s

//! [R] Molar Gas constant
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps 9.1e-7)
static double const gasConstantSI = 8.3144621; // J/mol/K

//! [k] BOLTZMANN'S CONSTANT == R/Na (JOULES/K)
// If this changes you msut update the Enumerated Temperature Type in UnitSystemEnusm.hh!
static double const boltzmannSI = 1.380648800E-23; // J K^-1

//! [e] ELECTRON CHARGE (COULOMBS)
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 2.2e-8)
// If this changes you msut update the Enumerated Temperature Type in UnitSystemEnusm.hh!
static double const electronChargeSI = 1.602176565e-19; // Amp / sec

//! [me] ELECTRON REST MASS (KG)	 s
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 4e-8)
static double const electronMassSI = 9.10938291e-31; // kg

//! [G] Gravitational Constant
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 1.2e-4)
static double const gravitationalConstantSI = 6.67384e-11; // m^3/kg/s^2

//! [g] Acceleration due to gravity (standard)
// Wikipedia (2013-12-3)
static double const accelerationFromGravitySI = 9.80665; // m/s^2

//! [mp] PROTON REST MASS (KG)
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 4.4e-8)
static double const protonMassSI = 1.672621777e-27; // kg

//---------------------------------------------------------------------------//
// DERIVED CONSTANTS
//  - constants derived from fundamental constants listed above
//  - constants with alternate units
//  - conversion factors
//---------------------------------------------------------------------------//

/*! \brief EV2K CONVERSION FACTOR FROM ELECTRON-VOLTS TO KELVIN (K/eV)
 *
 * (TEMPERATURE IN eV) * EV2K = (TEMPERATURE IN KELVIN)
 *
 * If this number is changed, you must also update the conversion factor found
 * in UniSystemUnums.hh.
 */
static double const EV2K = electronChargeSI / boltzmannSI;

/*! [sigma] STEFAN-BOLTZMANN CONSTANT (WATTS/(M**2-K**4)
 *
 * /f
 * \sigma_{SB} = \frac{2 \pi^5 k^4} {15 h^3 c^2}
 *             = 5.670373e-8
 * /f
 */
static double const stefanBoltzmannSI =
    static_cast<double>(2.0) * std::pow(PI, 5) * std::pow(boltzmannSI, 4) /
    (static_cast<double>(15.0) * std::pow(planckSI, 3) * std::pow(cLightSI, 2));

//! [F] Faraday constant == Na * e
static double const faradayConstantSI = AVOGADRO * electronChargeSI;

//! [mu0] Permeability of vacuum (free space)
static double const permeabilityOfVacuumSI = 4.0 * PI * 1.0e-7; // Henry/m

//! [epsi0] PERMITTIVITY OF FREE SPACE (F/M)
static double const permittivityOfFreeSpaceSI =
    1.0 / permeabilityOfVacuumSI / cLightSI / cLightSI; // Coloumb^2/J/m

//! [re] Classical electron radius (M)
static double const classicalElectronRadiusSI =
    std::pow(electronChargeSI, 2) / (4 * PI * permittivityOfFreeSpaceSI *
                                     electronMassSI * std::pow(cLightSI, 2));

} // end namespace rtt_units

#endif // __units_PhysicalConstantsSI_hh__

//---------------------------------------------------------------------------//
// end of units/PhysicalConstantsSI.hh
//---------------------------------------------------------------------------//
