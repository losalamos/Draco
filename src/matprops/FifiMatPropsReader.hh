//----------------------------------*-C++-*----------------------------------//
// FifiMatPropsReader.hh
// Randy M. Roberts
// Mon Apr 20 13:44:11 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_FifiMatPropsReader_hh__
#define __matprops_FifiMatPropsReader_hh__

#include "MaterialPropsReader.hh"
#include "FifiParser.hh"

#include "ds++/Mat.hh"

#include <string>
#include <vector>
#include <map>
#include <iosfwd>

namespace rtt_matprops {    

//===========================================================================//
// class FifiMatPropsReader - 
//
// Date created :
// Purpose      :  Concrete class for reading material property data
//                 from a Fifi file.  This data will be demanded of this
//                 object by a Material Props object.
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class FifiMatPropsReader : public MaterialPropsReader
{
    // NESTED CLASSES AND TYPEDEFS

  public:

    //========================================================================//
    // Nested struct:  MaterialDefinition
    //   This struct is used to further define the materials to be processed
    //   from the Fifi file.  The Fifi file has no means to determine the
    //   name of the material, nor the material's mass (in amus).
    //========================================================================//

    struct MaterialDefinition
    {
	std::string         name;
	MaterialId          matid;
	double              abar;
	MaterialDefinition() { }
	MaterialDefinition(const std::string &name_, MaterialId matid_,
			   double abar_)
	    : name(name_), matid(matid_), abar(abar_)
	{
	    // empty
	}
    };

  private:
    
    //========================================================================//
    // Nested struct:  MaterialInfo
    //   Each material will have its own temperature, density, and energy grid.
    //   Along with its mass in amus, these describe the material info.
    //========================================================================//

    struct MaterialInfo
    {
	std::string         name;
	MaterialId          matid;

	// Mean Atomic Weight (amu)
	double              abar; 

	std::vector<double> temperatureGrid;
	std::vector<double> densityGrid;
	std::vector<double> energyGrid;
	
	MaterialInfo()
	{
	    //*empty*
	}
	MaterialInfo(const std::string &name_, MaterialId matid_, double abar_)
	    : name(name_), matid(matid_), abar(abar_)
	{
	    //*empty*
	}
	const std::vector<double> &getTemperatureGrid() const
	{
	    return temperatureGrid;
	}
	const std::vector<double> &getDensityGrid() const
	{
	    return densityGrid;
	}
	const std::vector<double> &getEnergyGrid() const
	{
	    return energyGrid;
	}
	int getNumTemperatures() const { return temperatureGrid.size(); }
	int getNumDensities() const { return densityGrid.size(); }

	// Since the energy grid stores both upper and lower bounds
	// the number of entries in the energy grid is numGroups + 1.
	// (The upper bounds of group ig is the lower bounds of group ig+1.)

	int getNumGroups() const { return energyGrid.size() - 1; }
    };
    
    //========================================================================//
    // MatInfoMap:
    //   MatInfoMap is a map from the material id to the material info spec.
    //========================================================================//
    
    typedef std::map<MaterialId, MaterialInfo> MatInfoMap;

    // DATA

  private:

    // The FifiParser.
    // This object actuall reads the Fifi file, and is responsible
    // for locating data associated with materials and their keywords.

    FifiParser fifiParser;
    
    // This is the units converter from **most** of the units
    // used in the Fifi file to SI units.
    
    XTM::Units fileUnits;

    // This is the units converter from "file" units to
    // the user's output units.
    
    XTM::Units file2OutputUnits;

    // This is the actual mapping from material ids to the MaterialInfo spec.
    
    MatInfoMap materialInfoMap;

    // DISALLOWED DEFAULT METHODS
    
  private:

    // We won't allow anyone to copy the FifiMatPropsReader.
    
    FifiMatPropsReader(const FifiMatPropsReader &rhs);
    FifiMatPropsReader& operator=(const FifiMatPropsReader &rhs);

  public:

    // CREATORS
    
    //------------------------------------------------------------------------//
    // FifiMatPropsReader:
    //   This constructor takes a vector of MaterialDefinition's,
    //   the desired output Units, and the Fifi file's input stream.
    //   The vector of MaterialDefinition's decides which materials to "look at"
    //   from the file.  Any materials with ids not matching the "matid" field
    //   of any of matdefs' entries will not be processed, and therefore,
    //   ignored.  The "abar" field of the MaterialDefinition's are the mass
    //   in amus of the component material.
    //------------------------------------------------------------------------//

    FifiMatPropsReader(const std::vector<MaterialDefinition> &matdefs,
		       const XTM::Units &outputUnits_, std::istream &is_);

    // MANIPULATORS
    
    //------------------------------------------------------------------------//
    // getNextMaterial:
    //   This method doesn't really do much but verify that the material exists,
    //   and loads its name into the resulting string, and returns true.
    //   If the material does not exist it returns false.
    //------------------------------------------------------------------------//

    virtual bool getNextMaterial(MaterialId materialId_, std::string &name);

    //------------------------------------------------------------------------//
    // getTemperatureGrid:
    //   Return the temperature grid for the material given by materialId.
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //------------------------------------------------------------------------//

    virtual bool getTemperatureGrid(MaterialId materialId,
				    std::vector<double> &tempGrid_);
    
    //------------------------------------------------------------------------//
    // getDensityGrid:
    //   Return the density grid for the material given by materialId.
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //------------------------------------------------------------------------//

    virtual bool getDensityGrid(MaterialId materialId,
				std::vector<double> &densityGrid_);

    //------------------------------------------------------------------------//
    // getNumGroups:
    //   Return the number of energy groups for the material given by
    //   materialId.
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //------------------------------------------------------------------------//

    virtual bool getNumGroups(MaterialId materialId, int &numGroups);


    //------------------------------------------------------------------------//
    // getEnergyUpperbounds:
    //   Return the upper bounds of the specified group for the material
    //   given by materialId.
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //   Units: (temp)
    //------------------------------------------------------------------------//

    virtual bool getEnergyUpperbounds(MaterialId materialId, int group,
				      double &energyUpperbounds_);

    //------------------------------------------------------------------------//
    // getEnergyLowerbounds:
    //   Return the lower bounds of the specified group for the material
    //   given by materialId.
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //   Units: (temp)
    //------------------------------------------------------------------------//

    virtual bool getEnergyLowerbounds(MaterialId materialId, int group,
				      double &energyLowerbounds_);

    //------------------------------------------------------------------------//
    // getSigmaTotal:
    //   Return the total cross-section at the specified group for the
    //   material given by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //   If the fifi file parser cannot find the appropriate data keywords then
    //   return false.
    //   Units: (length^2/mass)
    //------------------------------------------------------------------------//

    virtual bool getSigmaTotal(MaterialId materialId, int group,
			       dsxx::Mat2<double> &data);

    //------------------------------------------------------------------------//
    // getSigmaAbsorption:
    //   Return the absorption cross-section at the specified group for the
    //   material given by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //   If the fifi file parser cannot find the appropriate data keywords then
    //   return false.
    //   Units: (length^2/mass)
    //------------------------------------------------------------------------//

    virtual bool getSigmaAbsorption(MaterialId materialId, int group,
				    dsxx::Mat2<double> &data);

    //------------------------------------------------------------------------//
    // getSigmaScattering:
    //   Return the scattering cross-section at the specified group for the
    //   material given by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //   If the fifi file parser cannot find the appropriate data keywords then
    //   return false.
    //   Units: (length^2/mass)
    //------------------------------------------------------------------------//

    virtual bool getSigmaScattering(MaterialId materialId, int group,
				    dsxx::Mat2<double> &data);

    //------------------------------------------------------------------------//
    // getSigmaEmission:
    //   Return the emission cross-section at the specified group for the
    //   material given by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //   If the fifi file parser cannot find the appropriate data keywords then
    //   return false.
    //   Units: (length^2/mass)
    //------------------------------------------------------------------------//

    virtual bool getSigmaEmission(MaterialId materialId, int group,
				  dsxx::Mat2<double> &data);

    //------------------------------------------------------------------------//
    // getElectronIonCoupling:
    //   Return the electron-ion coupling coefficient for the material given
    //   by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //   If the fifi file parser cannot find the appropriate data keywords then
    //   return false.
    //   Units: (energy/mass-time-temp)
    //------------------------------------------------------------------------//

    virtual bool getElectronIonCoupling(MaterialId materialId,
					dsxx::Mat2<double> &data);
	
    //------------------------------------------------------------------------//
    // getElectronConductionCoeff:
    //   Return the electron conduction coefficient for the material given
    //   by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //   If the fifi file parser cannot find the appropriate data keywords then
    //   return false.
    //   Units: (energy/length-time-temp-density)
    //------------------------------------------------------------------------//

    virtual bool getElectronConductionCoeff(MaterialId materialId,
					    dsxx::Mat2<double> &data);
	
    //------------------------------------------------------------------------//
    // getIonConductionCoeff:
    //   Return the ion conduction coefficient for the material given
    //   by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //   If the fifi file parser cannot find the appropriate data keywords then
    //   return false.
    //   Units: (energy/length-time-temp-density)
    //------------------------------------------------------------------------//

    virtual bool getIonConductionCoeff(MaterialId materialId,
				       dsxx::Mat2<double> &data);
	
    //------------------------------------------------------------------------//
    // getElectronSpecificHeat:
    //   Return the electron specific heat for the material given
    //   by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //   If the fifi file parser cannot find the appropriate data keywords then
    //   return false.
    //   Units: (energy/mass-temp)
    //------------------------------------------------------------------------//

    virtual bool getElectronSpecificHeat(MaterialId materialId,
					 dsxx::Mat2<double> &data);
	
    //------------------------------------------------------------------------//
    // getIonSpecificHeat:
    //   Return the ion specific heat for the material given
    //   by materialId.
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   The call to getMaterialInfo will throw an exception if the material
    //   is not found.
    //   If the fifi file parser cannot find the appropriate data keywords then
    //   return false.
    //   Units: (energy/mass-temp)
    //------------------------------------------------------------------------//

    virtual bool getIonSpecificHeat(MaterialId materialId,
				    dsxx::Mat2<double> &data);

    // ACCESSORS

    // *** none ***
    

    // IMPLEMENTATION

  private:

    //------------------------------------------------------------------------//
    // calcGridInfo:
    //    Go through all of the materials and fill in their temperature,
    //    density, and group energy grid information.
    //------------------------------------------------------------------------//

    void calcGridInfo();

    //------------------------------------------------------------------------//
    // hasMaterial:
    //    Returns true if material is found based on the mapping from
    //    the materialId.
    //------------------------------------------------------------------------//

    bool hasMaterial(MaterialId materialId) const;
    
    //------------------------------------------------------------------------//
    // getMaterialInfo:
    //   Return the MaterialInfo spec from the given materialId.
    //   Throw an exception if the material is not in the map.
    //------------------------------------------------------------------------//

    const MaterialInfo &getMaterialInfo(MaterialId materialId) const;

    //------------------------------------------------------------------------//
    // getSigma:
    //   Get a set of cross-section data, at the specified group for the
    //   material given by the MaterialInfo spec, from the Fifi file parser
    //   using the supplied keyword (e.g. "ramg", or "rsmg0").
    //   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
    //   Units: (length^2/mass)
    //------------------------------------------------------------------------//

    void getSigma(const MaterialInfo &matInfo, int group,
		  const std::string &keyword, dsxx::Mat2<double> &dataMat);

    //------------------------------------------------------------------------//
    // calcTemperatureDerivative:
    //    Calculate the temperature derivative of the nTemps x nDensities
    //    quantities.
    //------------------------------------------------------------------------//

    void calcTemperatureDerivative(MaterialId materialId,
				   const dsxx::Mat2<double> &data,
				   dsxx::Mat2<double> &derivative) const;
};

} // end of rtt_matprops namespace

#endif                          // __matprops_FifiMatPropsReader_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/FifiMatPropsReader.hh
//---------------------------------------------------------------------------//
