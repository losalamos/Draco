//----------------------------------*-C++-*----------------------------------//
// InterpedMaterialProps.cc
// Randy M. Roberts
// Wed Apr 15 08:44:49 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/InterpedMaterialProps.hh"

#include "matprops/BilinearInterpGrid.hh"
#include "matprops/BilinearInterpTable.hh"
#include "matprops/MaterialPropsReader.hh"

#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
using std::vector;

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM

typedef InterpedMaterialProps IMP;

IMP::InterpedMaterialProps(Units units_, const MaterialPropsReader &reader)
    : units(units_)
{
    int materialId;
    std::string materialName;
    
    while (reader.getNextMaterial(materialId, materialName))
    {
	// Make sure that this materialId
	// has not been used before.

	Assert(materials.find(materialId) == materials.end());

	vector<double> tempGrid;
	vector<double> densityGrid;
	reader.getTemperatureGrid(materialId, tempGrid);
	reader.getDensityGrid(materialId, densityGrid);

	SP<BilinearInterpGrid> spGrid(new BilinearInterpGrid(tempGrid,
							     densityGrid));
	int numGroups;
	reader.getNumGroups(materialId, numGroups);

	vector<double> energyUpperbounds(numGroups);
	vector<double> energyLowerbounds(numGroups);

	for (int i=0; i<numGroups; i++)
	{
	    reader.getEnergyUpperbounds(materialId, i, energyUpperbounds[i]);
	    reader.getEnergyLowerbounds(materialId, i, energyLowerbounds[i]);
	}

	dsxx::Mat2<double> data(tempGrid.size(), densityGrid.size(), 0.0);
	vector<BilinearInterpTable> tables(numGroups);

	for (int i=0; i<numGroups; i++)
	{
	    reader.getSigmaTotal(materialId, i, data);

	    tables[i] = BilinearInterpTable(spGrid, data);
	}

	GroupedTable sigmaTotal(energyUpperbounds, energyLowerbounds,
				tables);
	
	for (int i=0; i<numGroups; i++)
	{
	    reader.getSigmaAbsorption(materialId, i, data);

	    tables[i] = BilinearInterpTable(spGrid, data);
	}

	GroupedTable sigmaAbsorption(energyUpperbounds, energyLowerbounds,
				tables);
	
	for (int i=0; i<numGroups; i++)
	{
	    reader.getSigmaEmission(materialId, i, data);

	    tables[i] = BilinearInterpTable(spGrid, data);
	}

	GroupedTable sigmaEmission(energyUpperbounds, energyLowerbounds,
				tables);

	reader.getElectronIonCoupling(materialId, data);
	BilinearInterpTable electronIonCoupling(spGrid, data);
	
	reader.getElectronConductionCoeff(materialId, data);
	BilinearInterpTable electronConductionCoeff(spGrid, data);
	
	reader.getIonConductionCoeff(materialId, data);
	BilinearInterpTable ionConductionCoeff(spGrid, data);
	
	reader.getElectronSpecificHeat(materialId, data);
	BilinearInterpTable electronSpecificHeat(spGrid, data);
	
	reader.getIonSpecificHeat(materialId, data);
	BilinearInterpTable ionSpecificHeat(spGrid, data);

	typedef MatTabMap::value_type ValuePair;

	std::pair<MatTabMap::iterator, bool> retval;

	retval =
	    materials.insert(ValuePair(materialId,
				       MaterialTables(materialName,
						      spGrid,
						      sigmaTotal,
						      sigmaAbsorption,
						      sigmaEmission,
						      electronIonCoupling,
						      electronConductionCoeff,
						      ionConductionCoeff,
						      electronSpecificHeat,
						      ionSpecificHeat)));
	// Make sure the insertion succeeded.
	
	Assert(retval.second);
	
    }  // Loop over materials
}

template<class FT, class FT2>
IMP::MaterialStateField<FT> IMP::getMaterialState(const FT &density_,
						  const FT &electronTemp_,
						  const FT &ionTemp_,
						  const FT2 &matId_) const
{
    return MaterialStateField<FT>(*this, density_, electronTemp_,
				  ionTemp_, matId_);
}

template<class FT>
template<class FT2>
IMP::MaterialStateField<FT>::MaterialStateField(const IMP &matprops_,
						const FT &density_,
						const FT &electronTemp_,
						const FT &ionTemp_,
						const FT2 &matId_)
    : matprops(matprops_), theSize(density_.size()),
      density(density_.begin(), density_.end()),
      electronTemp(electronTemp_.begin(), electronTemp_.end()),
      ionTemp(ionTemp_.begin(), ionTemp_.end()),
      matId(matId_.begin(), matId_.end())
{
    Require(density_.size() == electronTemp_.size());
    Require(density_.size() == ionTemp_.size());
    Require(density_.size() == matId_.size());

    memento.reserve(size());
    
    for (int i=0; i < size(); i++)
    {
	const MaterialTables &mattabs = matprops.getMaterialTables(getMatId(i));
	
	const BilinearInterpGrid &grid = mattabs.getGrid();

	memento.push_back(grid.getMemento(getDensity(i), getElectronTemp(i)));
    }
}

//---------------------------------------------------------------------------//
// getValuesFromMatTable
//---------------------------------------------------------------------------//

template<class FT>
void IMP::getValuesFromMatTable(const MaterialStateField<FT> &matState,
				int group, PGroupedTable pTable,
				FT &results) const
{
    Require(&matState.matprops == this);
    Require(matState.size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < matState.size(); i++)
    {
	const MaterialTables &mattabs = getMaterialTables(matState.getMatId(i));

	const GroupedTable &groupTable = (mattabs.*pTable)();
	
	// groups are 1-based

	Require(group >= 1 && group <= groupTable.numGroups());

	const BilinearInterpTable &table = groupTable.getTable(group);

	table.interpolate(matState.getMemento(i), *resit++);
    }
}

//---------------------------------------------------------------------------//
// getValuesFromMatTable
//---------------------------------------------------------------------------//

template<class FT>
void IMP::getValuesFromMatTable(const MaterialStateField<FT> &matState,
				PBilinearInterpTable pTable, FT &results) const
{
    Require(&matState.matprops == this);
    Require(matState.size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < matState.size(); i++)
    {
	const MaterialTables &mattabs
	    = getMaterialTables(matState.getMatId(i));
	const BilinearInterpTable &table = (mattabs.*pTable)();
	
	table.interpolate(matState.getMemento(i), *resit++);
    }
}

template<class FT>
void IMP::getElectronTemperature(const MaterialStateField<FT> &matState,
				 FT &results) const
{
    Require(&matState.matprops == this);
    Require(matState.size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < matState.size(); i++)
	*resit++ = matState.getElectronTemp(i);
}

template<class FT>
void IMP::getIonTemperature(const MaterialStateField<FT> &matState,
			    FT &results) const
{
    Require(&matState.matprops == this);
    Require(matState.size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < matState.size(); i++)
	*resit++ = matState.getIonTemp(i);
}

template<class FT>
void IMP::getDensity(const MaterialStateField<FT> &matState, FT &results) const
{
    Require(&matState.matprops == this);
    Require(matState.size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < matState.size(); i++)
	*resit++ = matState.getDensity(i);
}

const IMP::MaterialTables &IMP::getMaterialTables(int matId) const
{
    MatTabMap::const_iterator matit = materials.find(matId);

    // Make sure the material exists.
	
    Assert(matit != materials.end());

    return (*matit).second;
}

IMP::MaterialTables &IMP::getMaterialTables(int matId)
{
    MatTabMap::iterator matit = materials.find(matId);

    // Make sure the material exists.
	
    Assert(matit != materials.end());

    return (*matit).second;
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of InterpedMaterialProps.cc
//---------------------------------------------------------------------------//
