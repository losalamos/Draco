//----------------------------------*-C++-*----------------------------------//
// InterpedMaterialProps.hh
// Randy M. Roberts
// Wed Apr 15 08:44:49 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __matprops_InterpedMaterialProps_hh__
#define __matprops_InterpedMaterialProps_hh__

#include "matprops/BilinearInterpGrid.hh"
#include "matprops/BilinearInterpTable.hh"
#include "ds++/SP.hh"
#include "3T/Units.hh"

#include <string>
#include <vector>
#include <map>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM

// Forward Reference

class MaterialPropsReader;


//===========================================================================//
// class InterpedMaterialProps - 
//
// Date created :
// Purpose      :
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class InterpedMaterialProps
{

    // NESTED CLASSES AND TYPEDEFS

  public:

    // Forward declaration
    
    template<class FT> class MaterialStateField;

  private:

    // Forward declaration
    
    struct GroupedTable;

    // Forward declaration
    
    struct MaterialTables;

    // Definine the type that contains a dictionary of material tables
    // keyed via the integer material id number.
    
    typedef std::map<int, MaterialTables> MatTabMap;


    // DATA

  private:
    
    Units units;

    // This is the dictionary of material tables
    // keyed via the integer material id number.
    
    MatTabMap materials;
    

    // CREATORS
    
  public:

    InterpedMaterialProps(Units units_, const MaterialPropsReader &reader);


    // MANIPULATORS
    
    // *** none ***


    // ACCESSORS

  public:

    Units getUnits() const { return units; }

    template<class FT, class FT2>
    MaterialStateField<FT> getMaterialState(const FT &density_,
					    const FT &electronTemp_,
					    const FT &ionTemp_,
					    const FT2 &matId_) const;    

    //------------------------------------------------------------------------//
    // getSigmaTotal:
    //------------------------------------------------------------------------//

    template<class FT>
    void getSigmaTotal(const MaterialStateField<FT> &matState,
		       int group, FT &results) const
    {
	getValuesFromMatTable(matState, group, &MaterialTables::sigmaTotal,
			      results);
    }

    template<class FT>
    void getSigmaTotal(const MaterialStateField<FT> &matState,
		       double group, FT &results) const
    {
	// Not yet implemented
	Assert(0);
    }

    template<class FT>
    void getSigmaAbsorption(const MaterialStateField<FT> &matState,
			    int group, FT &results) const
    {
	getValuesFromMatTable(matState, group,
			      &MaterialTables::getSigmaAbsorption, results);
    }

    template<class FT>
    void getSigmaAbsorption(const MaterialStateField<FT> &matState,
			    double group, FT &results) const
    {
	// Not yet implemented
	Assert(0);
    }

    template<class FT>
    void getSigmaEmission(const MaterialStateField<FT> &matState,
			  int group, FT &results) const
    {
	getValuesFromMatTable(matState, group,
			      &MaterialTables::getSigmaEmission, results);
    }

    template<class FT>
    void getSigmaEmission(const MaterialStateField<FT> &matState,
			  double group, FT &results) const
    {
	// Not yet implemented
	Assert(0);
    }

    template<class FT>
    void getElectronIonCoupling(const MaterialStateField<FT> &matState,
				FT &results) const
    {
	getValuesFromMatTable(matState, &MaterialTables::getElectronIonCoupling,
			      results);
    }

    template<class FT>
    void getElectronConductionCoeff(const MaterialStateField<FT> &matState,
				    FT &results) const
    {
	getValuesFromMatTable(matState,
			      &MaterialTables::getElectronConductionCoeff,
			      results);
    }

    template<class FT>
    void getIonConductionCoeff(const MaterialStateField<FT> &matState,
			       FT &results) const
    {
	getValuesFromMatTable(matState, &MaterialTables::getIonConductionCoeff,
			      results);
    }

    template<class FT>
    void getElectronSpecificHeat(const MaterialStateField<FT> &matState,
				 FT &results) const
    {
	getValuesFromMatTable(matState,
			      &MaterialTables::getElectronSpecificHeat,
			      results);
    }

    template<class FT>
    void getIonSpecificHeat(const MaterialStateField<FT> &matState,
			    FT &results) const
    {
	getValuesFromMatTable(matState, &MaterialTables::getIonSpecificHeat,
			      results);
    }
    
    template<class FT>
    void getElectronTemperature(const MaterialStateField<FT> &matState,
				FT &results) const;

    template<class FT>
    void getIonTemperature(const MaterialStateField<FT> &matState,
			   FT &results) const;

    template<class FT>
    void getDensity(const MaterialStateField<FT> &matState, FT &results) const;


    // IMPLEMENTATION

  private:
    
    const MaterialTables &getMaterialTables(int matId) const;
	
    MaterialTables &getMaterialTables(int matId);

    typedef GroupedTable (MaterialTables::*PGroupedTable)();
    
    template<class FT>
    void getValuesFromMatTable(const MaterialStateField<FT> &matState,
			       int group, PGroupedTable pTable,
			       FT &results) const;

    typedef BilinearInterpTable (MaterialTables::*PBilinearInterpTable)();
    
    template<class FT>
    void getValuesFromMatTable(const MaterialStateField<FT> &matState,
			       PBilinearInterpTable pTable,
			       FT &results) const;
};

//========================================================================//
// class InterpedMaterialProps::GroupedTable
//    This class is defined for material properties which vary
//    by group or continuous energy.
//    For energy, one must interpolate between groups (not yet implemented).
//========================================================================//

struct InterpedMaterialProps::GroupedTable
{
    // DATA

    int ngroups;
	
    // all of these vectors are ngroups long.
	
    std::vector<double>              energyUpperbounds;
    std::vector<double>              energyLowerbounds;
    std::vector<BilinearInterpTable> tables;

    // CREATORS

    GroupedTable() : ngroups(0) {};
    
    GroupedTable(const std::vector<double> &energyUpperbounds_,
		 const std::vector<double> &energyLowerbounds_,
		 const std::vector<BilinearInterpTable> tables_)
	: energyUpperbounds(energyUpperbounds_),
	  energyLowerbounds(energyLowerbounds_), tables(tables_)
    {
	ngroups = energyUpperbounds.size();
	Require(energyLowerbounds.size() == ngroups);
	Require(tables.size() == ngroups);
    }
    
    // ACCESSORS
	
    int numGroups() const { return ngroups; }

    const BilinearInterpTable &getTable(int groupNo) const
    {
	return tables[groupNo-1];
    }

    // MANIPULATORS
	
    BilinearInterpTable &getTable(int groupNo)
    {
	return tables[groupNo-1];
    }
};

//===========================================================================//
// class InterpedMaterialProps::MaterialTables
//===========================================================================//

struct InterpedMaterialProps::MaterialTables
{
    std::string              materialName;
    
    SP<BilinearInterpGrid>   spGrid;

    GroupedTable             sigmaTotal;
    GroupedTable             sigmaAbsorption;
    GroupedTable             sigmaEmission;
    BilinearInterpTable      electronIonCoupling;
    BilinearInterpTable      electronConductionCoeff;
    BilinearInterpTable      ionConductionCoeff;
    BilinearInterpTable      electronSpecificHeat;
    BilinearInterpTable      ionSpecificHeat;

    MaterialTables() {};
    
    MaterialTables(const std::string &materialName_,
		   const SP<BilinearInterpGrid> &spGrid_,
		   const GroupedTable &sigmaTotal_,
		   const GroupedTable &sigmaAbsorption_,
		   const GroupedTable &sigmaEmission_,
		   const BilinearInterpTable &electronIonCoupling_,
		   const BilinearInterpTable &electronConductionCoeff_,
		   const BilinearInterpTable &ionConductionCoeff_,
		   const BilinearInterpTable &electronSpecificHeat_,
		   const BilinearInterpTable &ionSpecificHeat_)
	: materialName(materialName_),
	  spGrid(spGrid_),
	  sigmaTotal(sigmaTotal_),
          sigmaAbsorption(sigmaAbsorption_),
	  sigmaEmission(sigmaEmission_),
          electronIonCoupling(electronIonCoupling_),
          electronConductionCoeff(electronConductionCoeff_),
          ionConductionCoeff(ionConductionCoeff_),
          electronSpecificHeat(electronSpecificHeat_),
          ionSpecificHeat(ionSpecificHeat_)
    {
	// **empty
    }

    const std::string &getMaterialName() const
    {
	return materialName;
    }
    
    const BilinearInterpGrid &getGrid() const
    {
	return *spGrid;
    }
    const GroupedTable &getSigmaTotal() const
    {
	return sigmaTotal;
    }
    const GroupedTable &getSigmaAbsorption() const
    {
	return sigmaAbsorption;
    }
    const GroupedTable &getSigmaEmission() const
    {
	return sigmaEmission;
    }
    const BilinearInterpTable &getElectronIonCoupling() const
    {
	return electronIonCoupling;
    }
    const BilinearInterpTable &getElectronConductionCoeff() const
    {
	return electronConductionCoeff;
    }
    const BilinearInterpTable &getIonConductionCoeff() const
    {
	return ionConductionCoeff;
    }
    const BilinearInterpTable &getElectronSpecificHeat() const
    {
	return electronSpecificHeat;
    }
    const BilinearInterpTable &getIonSpecificHeat() const
    {
	return ionSpecificHeat;
    }
};

//===========================================================================//
// class InterpedMaterialProps::MaterialStateField
//===========================================================================//

template<class FT>
class InterpedMaterialProps::MaterialStateField
{
    friend class InterpedMaterialProps;

    typedef BilinearInterpTable::Memento Memento;
    
  private:

    int theSize;
    
    const InterpedMaterialProps  &matprops;

    std::vector<Memento>         memento;
    
    std::vector<double>          density;
    std::vector<double>          electronTemp;
    std::vector<double>          ionTemp;
    std::vector<int>             matId;
	

  private:

    template<class FT2>
    MaterialStateField(const InterpedMaterialProps &matprops_,
		       const FT &density_, const FT &electronTemp_,
		       const FT &ionTemp_, const FT2 &matId_);

    int size() { return theSize; }

    const Memento &getMemento(int i) const { return memento[i]; }
    
    const double &getDensity(int i) const { return density[i]; }
    const double &getElectronTemp(int i) const { return electronTemp[i]; }
    const double &getIonTemp(int i) const { return ionTemp[i]; }
    const double &getMatId(int i) const { return matId[i]; }
};

END_NS_XTM  // namespace XTM

#endif                          // __matprops_InterpedMaterialProps_hh__

//---------------------------------------------------------------------------//
//                              end of matprops/InterpedMaterialProps.hh
//---------------------------------------------------------------------------//
