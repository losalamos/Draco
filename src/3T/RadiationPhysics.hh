//----------------------------------*-C++-*----------------------------------//
// RadiationPhysics.hh
// Randy M. Roberts
// Wed Mar 18 13:34:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_RadiationPhysics_hh__
#define __3T_RadiationPhysics_hh__

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
//===========================================================================//
// class RadiationPhysics - 
//
// Date created : Wed Mar 18 13:34:37 1998
// Purpose      : Calculate parameters and quantities needed
//                for radiation physics.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "3T/Units.hh"
#include "3T/PhysicalConstants.hh"

class RadiationPhysics
{

    // NESTED CLASSES AND TYPEDEFS
    // none

    // DATA

    Units units;  // The user's units.
    
  public:

    // CREATORS
    
    explicit RadiationPhysics(const Units &units_);

    // MANIPULATORS

    //------------------------------------------------------------------------//
    // setUnits:
    //     Allow user to change the units.
    //------------------------------------------------------------------------//

    void setUnits(const Units &units_)
    {
	units = units_;
    }

    // ACCESSORS

    //------------------------------------------------------------------------//
    // getUnits:
    //     Return the user's units to any interested party.
    //------------------------------------------------------------------------//

    const Units& getUnits() const { return units; }

    //------------------------------------------------------------------------//
    // getLightSpeed:
    //     Calculate the speed of light in the user's units
    //------------------------------------------------------------------------//

    inline double getLightSpeed() const;
    
    //------------------------------------------------------------------------//
    // getStefanBoltzmann:
    //     Calculate the Stefan-Boltzmann constant in the user's units.
    //------------------------------------------------------------------------//

    double getStefanBoltzmann() const;

    //------------------------------------------------------------------------//
    // getPlank:
    //   Calculate the planck function from the electron temperature.
    //   The units for input and output are specified by the units
    //   supplied to the class.
    //------------------------------------------------------------------------//

    template<class Field>
    void getPlanck(const Field &TElectron, Field &planckian) const;

    //------------------------------------------------------------------------//
    // getPlankTemperatureDerivative:
    //   Calculate the derivative of the planck function with respect to
    //   the electron temperature.
    //   The units for input and output are specified by the units
    //   supplied to the class.
    //------------------------------------------------------------------------//

    template<class Field>
    void getPlanckTemperatureDerivative(const Field &TElectron,
					Field &dPlanckdT) const;

  private:
    
    // IMPLEMENTATION
    // none
};

// INLINE DEFINITIONS

//------------------------------------------------------------------------//
// getLightSpeed:
//     Calculate the speed of light in the user's units
//------------------------------------------------------------------------//

inline double RadiationPhysics::getLightSpeed() const
{
    using XTM::PhysicalConstants::cLightSI;
    return units.SI2UserVelocity(cLightSI);
}

END_NS_XTM  // namespace XTM

#endif                          // __3T_RadiationPhysics_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/RadiationPhysics.hh
//---------------------------------------------------------------------------//
