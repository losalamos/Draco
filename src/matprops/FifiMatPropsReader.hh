//----------------------------------*-C++-*----------------------------------//
// FifiMatPropsReader.hh
// Randy M. Roberts
// Mon Apr 20 13:44:11 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_FifiMatPropsReader_hh__
#define __matprops_FifiMatPropsReader_hh__

#include "matprops/MaterialPropsReader.hh"
#include "matprops/FifiParser.hh"

#include "ds++/Mat.hh"

#include <string>
#include <vector>
#include <map>
#include <iosfwd>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    

//===========================================================================//
// class FifiMatPropsReader - 
//
// Date created :
// Purpose      :
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
    
    struct MaterialInfo
    {
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
	MaterialInfo(MaterialId matid_, double abar_)
	    : matid(matid_), abar(abar_)
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
	int getNumGroups() const { return energyGrid.size() - 1; }
    };
    
    typedef std::map<MaterialId, MaterialInfo> MatInfoMap;

    // DATA

  private:

    // The FifiParser

    FifiParser fifiParser;
    
    // This is the units converter from **most** of the units
    // used in the Fifi file to SI units.
    
    Units fileUnits;

    // This is the units converter from "file" units to
    // the user's output units.
    
    Units file2OutputUnits;

    MatInfoMap materialInfoMap;

    // DISALLOWED DEFAULT METHODS
    
  private:
    
    FifiMatPropsReader(const FifiMatPropsReader &rhs);
    FifiMatPropsReader& operator=(const FifiMatPropsReader &rhs);

  public:

    // CREATORS
    
    FifiMatPropsReader(const std::vector<MaterialDefinition> &matdefs,
		       const Units &outputUnits_, std::istream &is_);

    // MANIPULATORS
    
    virtual bool getNextMaterial(MaterialId materialId_, std::string &name);

    virtual bool getTemperatureGrid(MaterialId materialId,
				    std::vector<double> &tempGrid_);
    
    virtual bool getDensityGrid(MaterialId materialId,
				std::vector<double> &densityGrid_);

    virtual bool getNumGroups(MaterialId materialId, int &numGroups);


    virtual bool getEnergyUpperbounds(MaterialId materialId, int group,
				      double &energyUpperbounds_);

    virtual bool getEnergyLowerbounds(MaterialId materialId, int group,
				      double &energyLowerbounds_);

    virtual bool getSigmaTotal(MaterialId materialId, int group,
			       dsxx::Mat2<double> &data);

    virtual bool getSigmaAbsorption(MaterialId materialId, int group,
				    dsxx::Mat2<double> &data);

    virtual bool getSigmaEmission(MaterialId materialId, int group,
				  dsxx::Mat2<double> &data);

    virtual bool getElectronIonCoupling(MaterialId materialId,
					dsxx::Mat2<double> &data);
	
    virtual bool getElectronConductionCoeff(MaterialId materialId,
					    dsxx::Mat2<double> &data);
	
    virtual bool getIonConductionCoeff(MaterialId materialId,
				       dsxx::Mat2<double> &data);
	
    virtual bool getElectronSpecificHeat(MaterialId materialId,
					 dsxx::Mat2<double> &data);
	
    virtual bool getIonSpecificHeat(MaterialId materialId,
				    dsxx::Mat2<double> &data);

    // ACCESSORS

    // *** none **
    

    // IMPLEMENTATION

  private:

    void calcGridInfo();

    bool hasMaterial(MaterialId materialId) const;
    
    const MaterialInfo &getMaterialInfo(MaterialId materialId) const;

    void getSigma(const MaterialInfo &matInfo, int group,
		  const std::string &keyword, dsxx::Mat2<double> &dataMat);

    void calcTemperatureDerivative(MaterialId materialId,
				   const dsxx::Mat2<double> &data,
				   dsxx::Mat2<double> &derivative) const;
};

END_NS_XTM  // namespace XTM

#endif                          // __matprops_FifiMatPropsReader_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/FifiMatPropsReader.hh
//---------------------------------------------------------------------------//
