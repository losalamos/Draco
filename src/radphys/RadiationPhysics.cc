//----------------------------------*-C++-*----------------------------------//
// RadiationPhysics.cc
// Randy M. Roberts
// Wed Mar 18 13:34:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "radphys/RadiationPhysics.hh"

#include "radphys/ContainerTraits.hh"
#include "units/PhysicalConstants.hh"

#include <cmath>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//------------------------------------------------------------------------//
// RadiationPhysics:
//     Constructor based on Units.
//------------------------------------------------------------------------//

RadiationPhysics::RadiationPhysics(const Units &units_)
    : units(units_)
{
    // empty
}

//------------------------------------------------------------------------//
// getStefanBoltzmann:
//     Calculate the Stefan-Boltzmann constant in the user's units.
//------------------------------------------------------------------------//

double RadiationPhysics::getStefanBoltzmann() const
{
    using XTM::PhysicalConstants::stefanBoltzmannSI;
    using std::pow;

    // The Stefan-Boltzmann constant has the following SI units:
    //     W m^-2 K^-4 = kg s^-3 K^-4  (W - watts = kg m^2 s^-3)
    //
    // Therefore, we must convert temperature^-4, time^-3, and
    // mass^1 from SI units to the user's units, using the user's
    // Units.InvertXXX methods.
    
    const double stefanBoltzmann = units.InvertMass(
	units.ConvertTemperature(units.ConvertTime(stefanBoltzmannSI, 3), 4));

    return stefanBoltzmann;
}

//------------------------------------------------------------------------//
// getPlank:
//   Calculate the planck function from the electron temperature.
//   The units for input and output are specified by the units
//   supplied to the class.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getPlanck(const Field &TElectron,
				 Field &planckian) const
{
    Require(TElectron.size() == planckian.size());

    const double sb = getStefanBoltzmann();

    planckian = sb * TElectron*TElectron*TElectron*TElectron;
}

//------------------------------------------------------------------------//
// getPlankTemperatureDerivative:
//   Calculate the derivative of the planck function with respect to
//   the electron temperature.
//   The units for input and output are specified by the units
//   supplied to the class.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getPlanckTemperatureDerivative(const Field &TElectron,
						      Field &dplanckdT) const
{
    Require(TElectron.size() == dplanckdT.size());

    const double sb = getStefanBoltzmann();

    dplanckdT = 4.0 * sb * TElectron*TElectron*TElectron;
}

//------------------------------------------------------------------------//
// getElectIonCoupling:
//   Calculate the Electron-Ion Coupling function from the density,
//   electron temperature, free electrons per ion, average atomic weight.
//   The units for input and output are specified by the units
//   supplied to the class.
//   Exception: abar is in amu.
//------------------------------------------------------------------------//

template<class Field>
void RadiationPhysics::getElectIonCoupling(const Field &density,
					   const Field &TElectron,
					   const Field &z,
					   const Field &abar,
					   Field &electIonCoupling) const
{
    //.... CALCULATES ELECTRON-ION COUPLING COEFFICIENT
    //     USING DEVELOPMENT OF L.D. LANDAU (1936). THIS TREATMENT
    //     ASSUMES THAT THE MATERIAL IS AN IDEAL, CLASSICAL PLASMA, AND
    //     THAT THE ELECTRON AND ION TEMPERATURES ARE WITHIN A
    //     FACTOR OF 1836*abar OF EACH OTHER.

    //.... REFERENCE: LANDAU AND LIFSHITZ, "PHYSICAL KINETICS", 
    //                BUTTERWORTH-HEINEMANN, OXFORD 1995, SECTION 42

    //                MIYAMOTO, "PLASMA PHYSICS FOR NUCLEAR FUSION",
    //                MIT PRESS, 1989. SECTION 4.3

    Require(density.size() == electIonCoupling.size());
    Require(TElectron.size() == electIonCoupling.size());
    Require(z.size() == electIonCoupling.size());

    // eic_constant IS A FUNCTION OF VARIOUS PHYSICAL CONSTANTS IN SI UNITS:
    // eic_constant = (avog_num**3)*sqrt(mass_elect)*(elect_chg**4)*1.e9/ & 
    //                ( sqrt((2.*pi)**3)*(eps_zero**2)*sqrt(boltz_const) )

    // eic_constantSI: (amu**3-m**5-sqrt(K)/kg-s**3)
    
    const double eic_constantSI = 2.993811e+22;

    // Convert eic_constantSI to user's units using InvertXXX().

    const double eic_constant =
	units.InvertLength(
		units.InvertTemperature(
		    units.ConvertMass(units.ConvertTime(eic_constantSI, 3)),
		    0.5),
		5);

    // I must initialize the field with something.
    
    Field lambda_ei = TElectron;

    getElectIonCoulombLog(density, abar, z, TElectron, lambda_ei);
    
    electIonCoupling = eic_constant * lambda_ei *
	((z*z*z) / (abar*abar*abar)) * density / pow(TElectron, 1.5);
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of RadiationPhysics.cc
//---------------------------------------------------------------------------//
