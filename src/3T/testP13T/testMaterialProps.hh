//----------------------------------*-C++-*----------------------------------//
// testMaterialProps.hh
// Randy M. Roberts
// Tue Nov  3 15:46:13 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __3T_testP13T_testMaterialProps_hh__
#define __3T_testP13T_testMaterialProps_hh__

#include "matprops/MaterialProperties.hh"
#include "ds++/SP.hh"

namespace rtt_3T_testP13T
{

 using rtt_matprops::MaterialProperties;
 
 //===========================================================================//
 // class testMaterialProps - 
 //
 // Date created :
 // Purpose      :
 //
 // revision history:
 // -----------------
 // 0) original
 // 
 //===========================================================================//

 template<class MaterialStateFieldCC, class MaterialStateFieldFC, class MT>
 class testMaterialProps : public MaterialProperties<MT>
 {

     // NESTED CLASSES AND TYPEDEFS

     typedef typename MT::ccsf  ccsf;
     typedef typename MT::fcdsf fcdsf;

     // DATA

     const MaterialStateFieldCC &matstateCC;
     const MaterialStateFieldFC &matstateFC;
    
   public:

     // CREATORS
    
     testMaterialProps(const MaterialStateFieldCC &matstateCC_,
		       const MaterialStateFieldFC &matstateFC_)
	 : MaterialProperties<MT>(matstateCC_.getUnits()),
	   matstateCC(matstateCC_), matstateFC(matstateFC_)
     {
	 // empty
     }

     // MANIPULATORS
    
     // ACCESSORS

     void getElectronTemperature(ccsf &results) const
     {
	 matstateCC.getElectronTemperature(results);
     }
     
     void getElectronTemperature(fcdsf &results) const
     {
	 matstateFC.getElectronTemperature(results);
     }

     void getElectronConductionCoeff(ccsf &results) const
     {
	 matstateCC.getElectronConductionCoeff(results);
     }
     
     void getElectronConductionCoeff(fcdsf &results) const
     {
	 matstateFC.getElectronConductionCoeff(results);
     }

     void getElectronSpecificHeat(ccsf &results) const
     {
	 matstateCC.getElectronSpecificHeat(results);
     }
     
     void getElectronSpecificHeat(fcdsf &results) const
     {
	 matstateFC.getElectronSpecificHeat(results);
     }

     void getIonTemperature(ccsf &results) const
     {
	 matstateCC.getIonTemperature(results);
     }
     
     void getIonTemperature(fcdsf &results) const
     {
	 matstateFC.getIonTemperature(results);
     }

     void getIonConductionCoeff(ccsf &results) const
     {
	 matstateCC.getIonConductionCoeff(results);
     }
     
     void getIonConductionCoeff(fcdsf &results) const
     {
	 matstateFC.getIonConductionCoeff(results);
     }

     void getIonSpecificHeat(ccsf &results) const
     {
	 matstateCC.getIonSpecificHeat(results);
     }
     
     void getIonSpecificHeat(fcdsf &results) const
     {
	 matstateFC.getIonSpecificHeat(results);
     }

     void getSigmaTotal(int groupNo, ccsf &results) const
     {
	 matstateCC.getSigmaTotal(groupNo, results);
     }
     
     void getSigmaTotal(int groupNo, fcdsf &results) const
     {
	 matstateFC.getSigmaTotal(groupNo, results);
     }
     
     void getSigmaAbsorption(int groupNo, ccsf &results) const
     {
	 matstateCC.getSigmaAbsorption(groupNo, results);
     }

     void getSigmaAbsorption(int groupNo, fcdsf &results) const
     {
	 matstateFC.getSigmaAbsorption(groupNo, results);
     }
     
     void getSigmaScattering(int groupNo, ccsf &results) const
     {
	 matstateCC.getSigmaScattering(groupNo, results);
     }

     void getSigmaScattering(int groupNo, fcdsf &results) const
     {
	 matstateFC.getSigmaScattering(groupNo, results);
     }
     
     void getSigmaEmission(int groupNo, ccsf &results) const
     {
	 matstateCC.getSigmaEmission(groupNo, results);
     }
     
     void getSigmaEmission(int groupNo, fcdsf &results) const
     {
	 matstateFC.getSigmaEmission(groupNo, results);
     }

     void getElectronIonCoupling(ccsf &results) const
     {
	 matstateCC.getElectronIonCoupling(results);
     }
     
     void getElectronIonCoupling(fcdsf &results) const
     {
	 matstateFC.getElectronIonCoupling(results);
     }
     
   private:
    
     // IMPLEMENTATION
 };

} // end namespace rtt_3T_testP13T

#endif                          // __3T_testP13T_testMaterialProps_hh__

//---------------------------------------------------------------------------//
//                              end of 3T/testP13T/testMaterialProps.hh
//---------------------------------------------------------------------------//
