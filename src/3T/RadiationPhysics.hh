//----------------------------------*-C++-*----------------------------------//
// RadiationPhysics.hh
// Randy M. Roberts
// Wed Mar 18 13:34:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_RadiationPhysics_hh__
#define __3T_RadiationPhysics_hh__

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

    // DATA

    Units units;
    
  public:

    // CREATORS
    
    RadiationPhysics(const Units &units_);

    // MANIPULATORS

    setUnits(const Units &units_)
    {
	units = units_;
    }

    // ACCESSORS

    const Units& getUnits() const { return units; }

    inline double getLightSpeed() const;
    
    template<class Field>
    void getPlanck(const Field &TElectron, Field &planckian) const;

    template<class Field>
    void getPlanckTemperatureDerivative(const Field &TElectron,
					Field &dPlanckdT) const;

  private:
    
    // IMPLEMENTATION

};

// INLINE DEFINITIONS

inline double RadiationPhysics::getLightSpeed() const
{
    using PhysicalConstants::cLightSI;
    return units.SI2UserVelocity(cLightSI);
}

#endif                          // __3T_RadiationPhysics_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/RadiationPhysics.hh
//---------------------------------------------------------------------------//
