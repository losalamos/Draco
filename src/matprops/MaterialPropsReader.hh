//----------------------------------*-C++-*----------------------------------//
// MaterialPropsReader.hh
// Randy M. Roberts
// Mon Apr 20 10:27:42 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_MaterialPropsReader_hh__
#define __matprops_MaterialPropsReader_hh__

#include "3T/Units.hh"

#include "ds++/Mat.hh"

#include <string>
#include <vector>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

void funcYYY(dsxx::Mat2<double> &data);

BEGIN_NS_XTM

void funcXXX(dsxx::Mat2<double> &data);

//===========================================================================//
// class MaterialPropsReader - 
//    Abstract base class for reading material property data,
//    and giving that data to a Material Props object.
//===========================================================================//

class MaterialPropsReader
{
    // NESTED CLASSES AND TYPEDEFS

  public:

    typedef int MaterialId;

    // DATA

  private:
    
    Units outputUnits;
    
    // CREATORS
    
  public:

    MaterialPropsReader(const Units &outputUnits_)
	: outputUnits(outputUnits_)
    {
	// ** empty **
    }

    // MANIPULATORS
    
    virtual bool getNextMaterial(MaterialId materialId_, std::string &name) = 0;

    virtual bool getTemperatureGrid(MaterialId materialId,
				    std::vector<double> &tempGrid_) = 0;
    
    virtual bool getDensityGrid(MaterialId materialId,
				std::vector<double> &densityGrid_) = 0;

    virtual bool getNumGroups(MaterialId materialId, int &numGroups) = 0;


    virtual bool getEnergyUpperbounds(MaterialId materialId, int group,
				      double &energyUpperbounds_) = 0;

    virtual bool getEnergyLowerbounds(MaterialId materialId, int group,
				      double &energyLowerbounds_) = 0;

    virtual bool getSigmaTotal(MaterialId materialId, int group,
			       dsxx::Mat2<double> &data) = 0;

    virtual bool getSigmaAbsorption(MaterialId materialId, int group,
				    dsxx::Mat2<double> &data) = 0;

    virtual bool getSigmaEmission(MaterialId materialId, int group,
				  dsxx::Mat2<double> &data) = 0;

    virtual bool getElectronIonCoupling(MaterialId materialId,
					dsxx::Mat2<double> &data) = 0;
	
    virtual bool getElectronConductionCoeff(MaterialId materialId,
					    dsxx::Mat2<double> &data) = 0;
	
    virtual bool getIonConductionCoeff(MaterialId materialId,
				       dsxx::Mat2<double> &data) = 0;
	
    virtual bool getElectronSpecificHeat(MaterialId materialId,
					 dsxx::Mat2<double> &data) = 0;
	
    virtual bool getIonSpecificHeat(MaterialId materialId,
				    dsxx::Mat2<double> &data) = 0;

    // ACCESSORS

    const Units &getOutputUnits() const { return outputUnits; }

};

END_NS_XTM  // namespace XTM

#endif                          // __matprops_MaterialPropsReader_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/MaterialPropsReader.hh
//---------------------------------------------------------------------------//
