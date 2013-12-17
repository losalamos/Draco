//----------------------------------*-C++-*----------------------------------//
/*! \file   PhysicalConstants.hh
 *  \author Kelly Thompson
 *  \brief  Provide a single place where physical constants (pi, speed of
 *          light, etc) are defined.
 *  \date   Fri Nov 07 10:04:52 2003
 *  \note   Copyright (C) 2003-2013 Los Alamos National Security, LLC.
 *          All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __units_PhysicalConstants_hh__
#define __units_PhysicalConstants_hh__

#include "UnitSystem.hh"

//! \namespace rtt_units Namespace for units and physical constants
namespace rtt_units
{

// Mathematical constants
    
//! pi the ratio of a circle's circumference to its diameter (dimensionless)
static double const PI = 3.141592653589793238462643383279; 

// Euler's number (dimensionless)
static double const N_EULER = 2.7182818284590452353602874;

// Base physical constants in SI units:

//    m - meters, kg - kilograms, s - seconds, K - kelvin
//    W - watts, J - joules, C - coulombs, F - farads
//    mol - mole
    
//! [Na] AVOGADRO'S NUMBER ("entities"/MOLE)
// Wikipedia (2013-12-3) == NIST Codata 2010
static double const AVOGADRO          = 6.02214129e23;     // entities/mol
   
//! [h] Planck constant ( Energy seconds )
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 4.4e-8)
static double const planckSI          = 6.62606957e-34;    // J s

//! [R] Molar Gas constant
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps 9.1e-7)
static double const gasConstantSI     = 8.3144621;         // J/mol/K

//! [k] BOLTZMANN'S CONSTANT == R/Na (JOULES/K)
static double const boltzmannSI       = 1.380648800E-23;   // J K^-1

//! [e] ELECTRON CHARGE (COULOMBS)
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 2.2e-8)
static double const electronChargeSI  = 1.602176565e-19;    // Amp / sec

/*! \brief EV2K CONVERSION FACTOR FROM ELECTRON-VOLTS TO KELVIN (K/eV)
 *
 * (TEMPERATURE IN eV) * EV2K = (TEMPERATURE IN KELVIN)
 * 
 * If this number is changed, you must also update the conversion factor found
 * in UniSystemUnums.hh
 */
static double const EV2K              = electronChargeSI / boltzmannSI;

//! [c] SPEED OF LIGHT (M/S)
// exact value by NIST definition
static double const cLightSI          = 2.99792458e8;        // m s^-1

/*! [sigma] STEFAN-BOLTZMANN CONSTANT (WATTS/(M**2-K**4)
 *
 * /f
 * \sigma_{SB} = \frac{2 \pi^5 k^4} {15 h^3 c^2}
 *             = 5.670373e-8 
 * /f
 */
static double const stefanBoltzmannSI =
    2.0 * std::pow(PI,5) * std::pow(boltzmannSI,4)
    / (15.0 * std::pow(planckSI,3) * std::pow(cLightSI,2) );

//! [G] Gravitational Constant
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 1.2e-4)
static double const gravitationalConstantSI = 6.67384e-11; // m^3/kg/s^2

//! [g] Acceleration due to gravity (standard)
// Wikipedia (2013-12-3)
static double const accelerationFromGravitySI = 9.80665;     // m/s^2

//! [F] Faraday constant == Na * e
static double const faradayConstantSI = AVOGADRO * electronChargeSI;

//! [mu0] Permeability of vacuum (free space)
static double const permeabilityOfVacuumSI = 4.0 * PI * 1.0e-7; // Henry/m

//! [epsi0] PERMITTIVITY OF FREE SPACE (F/M)
static double const permittivityOfFreeSpaceSI =
    1.0 / permeabilityOfVacuumSI/cLightSI/cLightSI; // Coloumb^2/J/m

//! [me] ELECTRON REST MASS (KG)	 s
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 4e-8)
static double const electronMassSI    = 9.10938291e-31;        // kg

//! [mp] PROTON REST MASS (KG)
// Wikipedia (2013-12-3) == NIST Codata 2010 (eps = 4.4e-8)
static double const protonMassSI      = 1.672621777e-27;        // kg

//==============================================================================
/*!
 * \class PhysicalConstants
 *
 * \brief A class to define and encapsulate physical and mathematical
 * constants in the current UnitSystem.
 *
 * \sa rtt_units::UnitSystem
 * \sa rtt_units::UnitSystemType
 * \sa UnitSystemEnums.hh
 * \sa rtt_units::FundUnit
 * \sa Nuclide Chart
 *
 * Code Sample: 
 *
 * \verbatim
 * UnitSystem u( UnitSystem().getAstroUnits() );
 * PhyscialConstants pc( u );
 * double mass   = 1.0; // grams (AstroPhysics mass units).
 * double energy = mass * pc.c() * pc.c();  // E = mc^2 where energy is in
 *                                          // AstroPhysics units.
 * \endverbatim
 *
 * Example 
 * \example test/tstPhysicalConstants.cc
 * This is the unit regression test for the PhysicalConstants class.  It
 * demonstrates typical usage.
 */
//==============================================================================

class DLL_PUBLIC PhysicalConstants
{
  public:
    
    // Constructors.
    PhysicalConstants();
    explicit PhysicalConstants( UnitSystem const & u );

    //! \todo Make electronCharge and Avaragodo adjustable based on units.

    //! accesses PI (dimensionless)
    static double pi()               { return PI; };

    //! accesses conversion factor for eV to Kelvin (Temp/Energy)
    static double eV2K()             { return EV2K; } 

    //! accesses Avogadro's number (1/mole)
    static double avogadro()         { return AVOGADRO; }
    //! see avogadro()
    static double Na()               { return avogadro(); }

    //! access the Planck constant (units of energy-time)
    double planck()                  const { return d_planck; }
    //! see planck()
    double h()                       const { return planck(); }

    //! access the Gas constant (units of energy/mol/temp)
    double gasConstant()             const { return d_gasConstant; }
    //! see gasConstant()
    double R()                       const { return gasConstant(); }

    //! accesses the Boltzmann constant (Energy/Temp)
    double boltzmann()               const { return d_boltzmann; }
    //! see boltzmann()
    double k()                       const { return boltzmann(); }

    //! accesses the electron charge (Charge)
    double electronCharge()          const { return d_electronCharge; }
    //! see electronCharge()
    double e()                       const { return electronCharge(); }

    //! accesses the speed of light (units of velocity)
    double speedOfLight()            const { return d_cLight; }
    //! see speedOfLight()
    double c()                       const { return speedOfLight(); }

    //! accesses the StefanBoltzmann constant (Work/Area/Temp^4 )
    double stefanBoltzmann()         const { return d_stefanBoltzmann; }
    //! see stefanBoltzmann()
    double sigma()                   const { return stefanBoltzmann(); }

    //! accesses the gravitational constant
    double gravitationalConstant()   const { return d_gravitationalConstant; }
    //! see gravitationalConstant()
    double G()                       const { return gravitationalConstant(); }

    //! access the acceleration due to gravity (standard).
    double accelerationFromGravity() const { return d_accelerationFromGravity; }
    //! see accelerationFromGravity()
    double g()                       const { return accelerationFromGravity(); }

    //! access the Faraday constant
    double faradayConstant()         const { return d_faradayConstant; }
    //! see faradayConstant()
    double F()                       const { return faradayConstant(); }

    //! access the Permeability of vacuum (free space)
    double permeabilityOfVacuum()    const { return d_permeabilityOfVacuum; }
    //! see permeabilityOfVacuum()
    double mu0()                     const { return permeabilityOfVacuum(); }

    //! accesses the permittivity of free space (units of force/length)
    double permittivityOfFreeSpace() const { return d_permittivityOfFreeSpace; }
    //! see permittivityOfFreeSpace() 
    double epsi0()                   const { return permittivityOfFreeSpace(); }

    //! accesses the electron mass (units of mass)
    double electronMass()            const { return d_electronMass; }
    //! see electronMass()
    double me()                      const { return electronMass(); }

    //! accesses the proton mass (units of mass)
    double protonMass()              const { return d_protonMass; }
    //! see protonMass()
    double mp()                      const { return protonMass(); }

  private:

    // Base physical constants in SI units:

    //    m - meters, kg - kilograms, s - seconds, K - kelvin
    //    W - watts, J - joules, C - coulombs, F - farads
    //    mol - mole
    
    //! [h] Planck constant ( Energy seconds )
    double const d_planck;

    //! [R] Gas constant
    double const d_gasConstant;

    //! [k] BOLTZMANN'S CONSTANT == R/Na (JOULES/K)
    double const d_boltzmann;

    //! [e] ELECTRON CHARGE (COULOMBS)
    double const d_electronCharge;
    
    //! [c] SPEED OF LIGHT (M/S)
    double const d_cLight;

    //! [sigma] STEFAN-BOLTZMANN CONSTANT (WATTS/(M**2-K**4)
    double const d_stefanBoltzmann;
    
    //! [G] Gravitational Constant
    double const d_gravitationalConstant; // SI units are N / m^2 / kg^2 

    //! [g] Acceleration due to gravity (standard)
    double const d_accelerationFromGravity;

    //! [F] Faraday constant == Na * e
    double const d_faradayConstant; // SI units are coloumbs/mol

    //! [mu0] Permeability of vacuum (free space)
    double const d_permeabilityOfVacuum;

    //! [epsi0] PERMITTIVITY OF FREE SPACE (F/M)
    double const d_permittivityOfFreeSpace;

    //! [me] ELECTRON REST MASS (KG)	 
    double const d_electronMass;
    
    //! [mp] PROTON REST MASS (KG)
    double const d_protonMass;
    
}; // end class PhysicalConstants

} // end namespace rtt_units

#endif  // __units_PhysicalConstants_hh__

//---------------------------------------------------------------------------//
// end of units/PhysicalConstants.hh
//---------------------------------------------------------------------------//
