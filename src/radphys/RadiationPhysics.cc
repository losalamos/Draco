//----------------------------------*-C++-*----------------------------------//
// RadiationPhysics.cc
// Randy M. Roberts
// Wed Mar 18 13:34:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "radphys/RadiationPhysics.hh"

#include "traits/ContainerTraits.hh"
#include "units/PhysicalConstants.hh"

#include <cmath>
#include <limits>

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

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of RadiationPhysics.cc
//---------------------------------------------------------------------------//
