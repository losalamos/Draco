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
#include <vector>
#include <map>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM
    
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

    template<class FT> class MaterialStateField;

  private:

    struct GroupedTable
    {
	std::vector<double>              energyUpperbounds;
	std::vector<double>              energyLowerbounds;
	std::vector<BilinearInterpTable> tables;
	
	int numGroups() const { return tables.size(); }

	const BilinearInterpTable &getTable(int groupNo) const
	{
	    return tables[groupNo-1];
	}
	BilinearInterpTable &getTable(int groupNo)
	{
	    return tables[groupNo-1];
	}
    };
	
    struct MaterialTables
    {
	SP<BilinearInterpGrid>   spGrid;

	GroupedTable             sigmaTotal;
	GroupedTable             sigmaAbsorption;
	GroupedTable             sigmaEmission;
	BilinearInterpTable      electronIonCoupling;
	BilinearInterpTable      electronConductionCoeff;
	BilinearInterpTable      ionConductionCoeff;
	BilinearInterpTable      electronSpecificHeat;
	BilinearInterpTable      ionSpecificHeat;

	const BilinearInterpGrid &getGrid() const { return *spGrid; }
    };

    typedef std::map<int, MaterialTables> MatTabMap;

    // DATA

  private:
    
    Units units;

    MatTabMap materials;
    
  public:

    // CREATORS
    
    InterpedMaterialProps(Units units_, MaterialPropsFactory &factory);

    // MANIPULATORS
    
    // *** none ***

    // ACCESSORS

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
	getValuesFromMatTable(matState, group, &MaterialTables::sigmaAbsorption,
			      results);
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
	getValuesFromMatTable(matState, group, &MaterialTables::sigmaEmission,
			      results);
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
	getValuesFromMatTable(matState, &MaterialTables::electronIonCoupling,
			      results);
    }

    template<class FT>
    void getElectronConductionCoeff(const MaterialStateField<FT> &matState,
				    FT &results) const
    {
	getValuesFromMatTable(matState,
			      &MaterialTables::electronConductionCoeff,
			      results);
    }

    template<class FT>
    void getIonConductionCoeff(const MaterialStateField<FT> &matState,
			       FT &results) const
    {
	getValuesFromMatTable(matState, &MaterialTables::ionConductionCoeff,
			      results);
    }

    template<class FT>
    void getElectronSpecificHeat(const MaterialStateField<FT> &matState,
				 FT &results) const
    {
	getValuesFromMatTable(matState, &MaterialTables::electronSpecificHeat,
			      results);
    }

    template<class FT>
    void getIonSpecificHeat(const MaterialStateField<FT> &matState,
			    FT &results) const
    {
	getValuesFromMatTable(matState, &MaterialTables::ionSpecificHeat,
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

  private:
    
    // IMPLEMENTATION

    const MaterialTables &getMaterialTables(int matId) const
    {
	MatTabMap::const_iterator matit = materials.find(matId);

	// Make sure the material exists.
	
	Assert(matit != materials.end());

	return (*matit).second;
    }
	
    MaterialTables &getMaterialTables(int matId)
    {
	MatTabMap::iterator matit = materials.find(matId);

	// Make sure the material exists.
	
	Assert(matit != materials.end());

	return (*matit).second;
    }

    typedef MaterialTables::*GroupedTable PGroupedTable;
    
    template<class FT>
    void getValuesFromMatTable(const MaterialStateField<FT> &matState,
			       int group, PGroupedTable pTable,
			       FT &results) const;

    typedef MaterialTables::*BilinearInterpTable PBilinearInterpTable;
    
    template<class FT>
    void getValuesFromMatTable(const MaterialStateField<FT> &matState,
			       PBilinearInterpTable pTable,
			       FT &results) const;
};

//===========================================================================//
// class MaterialStateField
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
