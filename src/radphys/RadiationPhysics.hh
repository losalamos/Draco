//----------------------------------*-C++-*----------------------------------//
// RadiationPhysics.hh
// Randy M. Roberts
// Wed Mar 18 13:34:37 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __radphys_RadiationPhysics_hh__
#define __radphys_RadiationPhysics_hh__

#include "units/Units.hh"
#include "units/PhysicalConstants.hh"

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
    // getRadConstant:
    //     Calculate the constant, a, for Planck = (a c / 4pi) T^4,
    //     in the user's units.
    //     a = 4 stefanBoltzmann / c
    //------------------------------------------------------------------------//

    inline double getRadConstant() const;

    //------------------------------------------------------------------------//
    // getPlanck:
    //   Calculate the planck function from the electron temperature.
    //   The units for input and output are specified by the units
    //   supplied to the class.
    //   planck = stefanBoltzmann / pi * T^4
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

    //------------------------------------------------------------------------//
    // getElectIonCoupling:
    //   Calculate the Electron-Ion Coupling function from the density,
    //   electron temperature, free electrons per ion, average atomic weight.
    //   The units for input and output are specified by the units
    //   supplied to the class.
    //   Exception: abar is in amu.
    //------------------------------------------------------------------------//

    template<class Field>
    void getElectIonCoupling(const Field &density, const Field &TElectron,
			     const Field &z, double abar,
			     Field &electIonCoupling) const;

    //------------------------------------------------------------------------//
    // getElectIonCoulombLog:
    //   Calculate the Electron-Ion Coulomb Log function from the density,
    //   electron temperature, free electrons per ion, average atomic weight.
    //   The units for input and output are specified by the units
    //   supplied to the class.
    //   Exception: abar is in amu.
    //------------------------------------------------------------------------//

    template<class Field>
    void getElectIonCoulombLog(const Field &density, const Field &TElectron,
			     const Field &z, double abar,
			     Field &lambda_ei) const;

    //------------------------------------------------------------------------//
    // getElectronConductionCoeff:
    //   Calculate the Electron Conduction Coefficient from the density,
    //   electron temperature, free electrons per ion, average atomic weight.
    //   The units for input and output are specified by the units
    //   supplied to the class.
    //   Exception: abar is in amu.
    //------------------------------------------------------------------------//

    template<class Field>
    void getElectronConductionCoeff(const Field &density,
				    const Field &TElectron,
				    const Field &z, double abar,
				    Field &electCondCoeff) const;

    //------------------------------------------------------------------------//
    // getElectElectCoulombLog:
    //   Calculate the Electron-Electron Coulomb Log function from the density,
    //   electron temperature, free electrons per ion, average atomic weight.
    //   The units for input and output are specified by the units
    //   supplied to the class.
    //   Exception: abar is in amu.
    //------------------------------------------------------------------------//

    template<class Field>
    void getElectElectCoulombLog(const Field &density,
				 const Field &TElectron,
				 const Field &z, double abar,
				 Field &lambda_ee) const;

    //------------------------------------------------------------------------//
    // getElectronGamma0:
    //    Returns a dimensionless constant for use in
    //    calculation of electron thermal conduction
    //    coefficeint. Correlation developed by
    //    C. Cranfill, LANL.
    //------------------------------------------------------------------------//

    template<class Field>
    void getElectronGamma0(const Field &z, Field &gamma0) const;

    //------------------------------------------------------------------------//
    // getIonConductionCoeff:
    //   Calculate the Ion Conduction Coefficient from the density,
    //   electron temperature, free electrons per ion, average atomic weight.
    //   The units for input and output are specified by the units
    //   supplied to the class.
    //   Exception: abar is in amu.
    //------------------------------------------------------------------------//

    template<class Field>
    void getIonConductionCoeff(const Field &density,
			       const Field &TElectron,
			       const Field &z, double abar,
			       Field &ionCondCoeff) const;

    //------------------------------------------------------------------------//
    // getIonIonCoulombLog:
    //   Calculate the Ion-Ion Coulomb Log function from the density,
    //   electron temperature, free electrons per ion, average atomic weight.
    //   The units for input and output are specified by the units
    //   supplied to the class.
    //   Exception: abar is in amu.
    //------------------------------------------------------------------------//

    template<class Field>
    void getIonIonCoulombLog(const Field &density,
				 const Field &TElectron,
				 const Field &z, double abar,
				 Field &lambda_ee) const;

    //------------------------------------------------------------------------//
    // getIonGamma0:
    //    Returns a dimensionless constant for use in
    //    calculation of ion thermal conduction
    //    coefficeint. Correlation developed by
    //    C. Cranfill, LANL.
    //------------------------------------------------------------------------//

    template<class Field>
    void getIonGamma0(const Field &z, const Field &TElect,
		      const Field &Tion, double abar, Field &gamma0) const;

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

    // Convert from SI to user's units with the InvertXXX method.
    
    return units.InvertVelocity(cLightSI);
}

//------------------------------------------------------------------------//
// getRadConstant:
//     Calculate the constant, a, for Planck = (a c / 4pi) T^4,
//     in the user's units.
//     a = 4 stefanBoltzmann / c
//------------------------------------------------------------------------//

inline double RadiationPhysics::getRadConstant() const
{
    return 4.0 / getLightSpeed() * getStefanBoltzmann();
}

END_NS_XTM  // namespace XTM

#endif                          // __radphys_RadiationPhysics_hh__

//---------------------------------------------------------------------------//
//                              end of radphys/RadiationPhysics.hh
//---------------------------------------------------------------------------//
