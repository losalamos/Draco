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

IMP::InterpedMaterialProps(Units units_, MaterialPropsFactor &factory)
    : units(units_)
{
    vector<double> energyGrid = factory.getEnergyGrid();
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

template<class FT, class FT2>
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

	const GroupedTable &groupTable = mattabs.*pTable;
	
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
	const BilinearInterpTable &table = mattabs.*pTable;
	
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

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of InterpedMaterialProps.cc
//---------------------------------------------------------------------------//
