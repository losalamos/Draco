//----------------------------------*-C++-*----------------------------------//
// MaterialPropsReader.hh
// Randy M. Roberts
// Mon Apr 20 10:27:42 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_MaterialPropsReader_hh__
#define __matprops_MaterialPropsReader_hh__

#include "units/Units.hh"

#include "ds++/Mat.hh"

#include <string>
#include <vector>

void funcYYY(rtt_dsxx::Mat2<double> &data);

namespace rtt_matprops {

void funcXXX(rtt_dsxx::Mat2<double> &data);

//===========================================================================//
// class MaterialPropsReader - 
//    Abstract base class for reading material property data,
//    and giving that data to a Material Props object.
//===========================================================================//

class MaterialPropsReader
{
    // NESTED CLASSES AND TYPEDEFS

  public:

    // How to reference a material.
    
    typedef int MaterialId;

    // DATA

  private:

    // The desired units to convert all of the data.
    
    rtt_units::Units outputUnits;
    
    // CREATORS
    
  public:

    MaterialPropsReader(const rtt_units::Units &outputUnits_)
	: outputUnits(outputUnits_)
    {
	// ** empty **
    }

    // MANIPULATORS
    
    //------------------------------------------------------------------------//
    // getNextMaterial:
    //   This method doesn't really do much but verify that the material exists,
    //   and loads its name into the resulting string, and returns true.
    //   If the material does not exist it returns false.
    //------------------------------------------------------------------------//

    virtual bool getNextMaterial(MaterialId materialId_, std::string &name) = 0;

    //------------------------------------------------------------------------//
    // getTemperatureGrid:
    //   Return the temperature grid for the material given by materialId.
    //------------------------------------------------------------------------//

    virtual bool getTemperatureGrid(MaterialId materialId,
				    std::vector<double> &tempGrid_) = 0;
    
    //------------------------------------------------------------------------//
    // getDensityGrid:
    //   Return the density grid for the material given by materialId.
    //------------------------------------------------------------------------//

    virtual bool getDensityGrid(MaterialId materialId,
				std::vector<double> &densityGrid_) = 0;

    //------------------------------------------------------------------------//
    // getNumGroups:
    //   Return the number of energy groups for the material given by
    //   materialId.
    //------------------------------------------------------------------------//

    virtual bool getNumGroups(MaterialId materialId, int &numGroups) = 0;


    //------------------------------------------------------------------------//
    // getEnergyUpperbounds:
    //   Return the upper bounds of the specified group for the material
    //   given by materialId.
    //   Units: (temp)
    //------------------------------------------------------------------------//

    virtual bool getEnergyUpperbounds(MaterialId materialId, int group,
				      double &energyUpperbounds_) = 0;

    //------------------------------------------------------------------------//
    // getEnergyLowerbounds:
    //   Return the lower bounds of the specified group for the material
    //   given by materialId.
    //   Units: (temp)
    //------------------------------------------------------------------------//

    virtual bool getEnergyLowerbounds(MaterialId materialId, int group,
				      double &energyLowerbounds_) = 0;

    //------------------------------------------------------------------------//
    // getSigmaTotal:
    //   Return the total cross-section at the specified group for the
    //   material given by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   If the implementation cannot find the appropriate data then
    //   return false.
    //   Units: (length^2/mass)
    //------------------------------------------------------------------------//

    virtual bool getSigmaTotal(MaterialId materialId, int group,
			       rtt_dsxx::Mat2<double> &data) = 0;

    //------------------------------------------------------------------------//
    // getSigmaAbsorption:
    //   Return the absorption cross-section at the specified group for the
    //   material given by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   If the implementation cannot find the appropriate data then
    //   return false.
    //   Units: (length^2/mass)
    //------------------------------------------------------------------------//

    virtual bool getSigmaAbsorption(MaterialId materialId, int group,
				    rtt_dsxx::Mat2<double> &data) = 0;

    //------------------------------------------------------------------------//
    // getSigmaScattering:
    //   Return the scattering cross-section at the specified group for the
    //   material given by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   If the implementation cannot find the appropriate data then
    //   return false.
    //   Units: (length^2/mass)
    //------------------------------------------------------------------------//

    virtual bool getSigmaScattering(MaterialId materialId, int group,
				    rtt_dsxx::Mat2<double> &data) = 0;

    //------------------------------------------------------------------------//
    // getSigmaEmission:
    //   Return the emission cross-section at the specified group for the
    //   material given by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   If the implementation cannot find the appropriate data then
    //   return false.
    //   Units: (length^2/mass)
    //------------------------------------------------------------------------//

    virtual bool getSigmaEmission(MaterialId materialId, int group,
				  rtt_dsxx::Mat2<double> &data) = 0;

    //------------------------------------------------------------------------//
    // getElectronIonCoupling:
    //   Return the electron-ion coupling coefficient for the material given
    //   by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   If the implementation cannot find the appropriate data then
    //   return false.
    //   Units: (energy/mass-time-temp)
    //------------------------------------------------------------------------//

    virtual bool getElectronIonCoupling(MaterialId materialId,
					rtt_dsxx::Mat2<double> &data) = 0;
	
    //------------------------------------------------------------------------//
    // getElectronConductionCoeff:
    //   Return the electron conduction coefficient for the material given
    //   by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   If the implementation cannot find the appropriate data then
    //   return false.
    //   Units: (energy/length-time-temp-density)
    //------------------------------------------------------------------------//

    virtual bool getElectronConductionCoeff(MaterialId materialId,
					    rtt_dsxx::Mat2<double> &data) = 0;
	
    //------------------------------------------------------------------------//
    // getIonConductionCoeff:
    //   Return the ion conduction coefficient for the material given
    //   by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   If the implementation cannot find the appropriate data then
    //   return false.
    //   Units: (energy/length-time-temp-density)
    //------------------------------------------------------------------------//

    virtual bool getIonConductionCoeff(MaterialId materialId,
				       rtt_dsxx::Mat2<double> &data) = 0;
	
    //------------------------------------------------------------------------//
    // getElectronSpecificHeat:
    //   Return the electron specific heat for the material given
    //   by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   If the implementation cannot find the appropriate data then
    //   return false.
    //   Units: (energy/mass-temp)
    //------------------------------------------------------------------------//

    virtual bool getElectronSpecificHeat(MaterialId materialId,
					 rtt_dsxx::Mat2<double> &data) = 0;
	
    //------------------------------------------------------------------------//
    // getIonSpecificHeat:
    //   Return the ion specific heat for the material given
    //   by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   If the implementation cannot find the appropriate data then
    //   return false.
    //   Units: (energy/mass-temp)
    //------------------------------------------------------------------------//

    virtual bool getIonSpecificHeat(MaterialId materialId,
				    rtt_dsxx::Mat2<double> &data) = 0;

    // ACCESSORS

    //------------------------------------------------------------------------//
    // getOutputUnits:
    //   Return the desired units for output.
    //------------------------------------------------------------------------//

    const rtt_units::Units &getOutputUnits() const { return outputUnits; }

};

} // end of rtt_matprops namespace

#endif                          // __matprops_MaterialPropsReader_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/MaterialPropsReader.hh
//---------------------------------------------------------------------------//
