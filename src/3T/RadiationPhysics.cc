//----------------------------------*-C++-*----------------------------------//
// RadiationPhysics.cc
// Randy M. Roberts
// Wed Mar 18 13:34:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "3T/RadiationPhysics.hh"
#include "3T/PhysicalConstants.hh"
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
    // mass^1 from SI units to the user's units.
    
    const double stefanBoltzmann =
	stefanBoltzmannSI * pow(units.getTemperatureConversion(), 4)
	* pow(units.getTimeConversion(), 3) / units.getMassConversion();

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

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of RadiationPhysics.cc
//---------------------------------------------------------------------------//
