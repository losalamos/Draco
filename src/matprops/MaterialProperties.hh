//----------------------------------*-C++-*----------------------------------//
// MaterialProperties.hh
// Randy M. Roberts
// Tue Nov  3 15:01:42 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_MaterialProperties_hh__
#define __matprops_MaterialProperties_hh__

#include "units/Units.hh"

namespace rtt_matprops
{
 //===========================================================================//
 // class MaterialProperties - 
 //
 // Date created :
 // Purpose      :
 //
 // revision history:
 // -----------------
 // 0) original
 // 
 //===========================================================================//

 template<class MT>
 class MaterialProperties
 {

     // NESTED CLASSES AND TYPEDEFS

     typedef MT MeshType;
     typedef typename MT::ccsf ccsf;
     typedef typename MT::fcdsf fcdsf;

     // DATA
    
     XTM::Units units;

   public:

     // CREATORS
    
     MaterialProperties(const XTM::Units &units_)
	 : units(units_)
     {
	 // empty
     }

     // MANIPULATORS
    
     // ACCESSORS

     const XTM::Units &getUnits() const { return units; }

     virtual void getElectronTemperature(ccsf &results) const = 0;
     virtual void getElectronTemperature(fcdsf &results) const = 0;

     virtual void getElectronConductionCoeff(ccsf &results) const = 0;
     virtual void getElectronConductionCoeff(fcdsf &results) const = 0;

     virtual void getElectronSpecificHeat(ccsf &results) const = 0;
     virtual void getElectronSpecificHeat(fcdsf &results) const = 0;

     virtual void getIonTemperature(ccsf &results) const = 0;
     virtual void getIonTemperature(fcdsf &results) const = 0;

     virtual void getIonConductionCoeff(ccsf &results) const = 0;
     virtual void getIonConductionCoeff(fcdsf &results) const = 0;

     virtual void getIonSpecificHeat(ccsf &results) const = 0;
     virtual void getIonSpecificHeat(fcdsf &results) const = 0;

     virtual void getSigmaTotal(int groupNo, ccsf &results) const = 0;
     virtual void getSigmaTotal(int groupNo, fcdsf &results) const = 0;
     
     virtual void getSigmaAbsorption(int groupNo, ccsf &results) const = 0;
     virtual void getSigmaAbsorption(int groupNo, fcdsf &results) const = 0;
     
     virtual void getSigmaEmission(int groupNo, ccsf &results) const = 0;
     virtual void getSigmaEmission(int groupNo, fcdsf &results) const = 0;

     virtual void getElectronIonCoupling(ccsf &results) const = 0;
     virtual void getElectronIonCoupling(fcdsf &results) const = 0;
     
   private:
    
     // IMPLEMENTATION
 };

} // end namespace rtt_matprops

#endif                          // __matprops_MaterialProperties_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/MaterialProperties.hh
//---------------------------------------------------------------------------//
