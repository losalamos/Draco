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

#include <string>
using std::string;
#include <vector>
using std::vector;

#include <stdexcept>
#include <strstream>

#ifndef BEGIN_NS_XTM
#define BEGIN_NS_XTM namespace XTM  {
#define END_NS_XTM }
#endif

BEGIN_NS_XTM

typedef InterpedMaterialProps IMP;

//===========================================================================//
// InterpedMaterialProps Methods
//===========================================================================//

//---------------------------------------------------------------------------//
// InterpedMaterialProps Constructor:
//    Create an InterpedMaterialProps from a units object and a concrete
//    object derived from MaterialPropsReader.
//---------------------------------------------------------------------------//

IMP::InterpedMaterialProps(const vector<int> &materialIds,
			   MaterialPropsReader &reader)
    : units(reader.getOutputUnits())
{
    
    for (int i=0; i<materialIds.size(); i++)
    {
	std::string materialName;
	int materialId = materialIds[i];

	bool found = reader.getNextMaterial(materialId, materialName);

	if (!found)
	{
	    std::ostrstream os;
	    os << "InterpedMaterialProps ctor: Unable to find material "
	       << materialId << " in reader." << std::ends;
	    throw std::runtime_error(os.str());
	}
	
	// Make sure that this materialId
	// has not been used before.

	if (materials.find(materialId) != materials.end())
	{
	    throw std::runtime_error(string("InterpedMaterialProps ctor: ")
				     + string("Duplicate material ids."));
	}

	vector<double> tempGrid;
	vector<double> densityGrid;
	reader.getTemperatureGrid(materialId, tempGrid);
	reader.getDensityGrid(materialId, densityGrid);

	bool errorOnOutOfBounds = false;

	SP<BilinearInterpGrid>
	    spGrid(new BilinearInterpGrid(tempGrid, densityGrid,
					  errorOnOutOfBounds));
	
	int numGroups;
	reader.getNumGroups(materialId, numGroups);

	vector<double> energyUpperbounds(numGroups);
	vector<double> energyLowerbounds(numGroups);

	for (int i=1; i<=numGroups; i++)
	{
	    reader.getEnergyUpperbounds(materialId, i, energyUpperbounds[i-1]);
	    reader.getEnergyLowerbounds(materialId, i, energyLowerbounds[i-1]);
	}

	dsxx::Mat2<double> data(tempGrid.size(), densityGrid.size(), 0.0);
	vector<BilinearInterpTable> tables(numGroups);

	for (int i=1; i<=numGroups; i++)
	{
	    if (reader.getSigmaTotal(materialId, i, data))
		tables[i-1] = BilinearInterpTable(spGrid, data);
	}

	GroupedTable sigmaTotal(energyUpperbounds, energyLowerbounds,
				tables);
	
	for (int i=1; i<=numGroups; i++)
	{
	    if (reader.getSigmaAbsorption(materialId, i, data))
		tables[i-1] = BilinearInterpTable(spGrid, data);
	}

	GroupedTable sigmaAbsorption(energyUpperbounds, energyLowerbounds,
				tables);
	
	for (int i=1; i<=numGroups; i++)
	{
	    if (reader.getSigmaEmission(materialId, i, data))
		tables[i-1] = BilinearInterpTable(spGrid, data);
	}

	GroupedTable sigmaEmission(energyUpperbounds, energyLowerbounds,
				   tables);

	BilinearInterpTable electronIonCoupling;
	if (reader.getElectronIonCoupling(materialId, data))
	    electronIonCoupling = BilinearInterpTable(spGrid, data);
	
	BilinearInterpTable electronConductionCoeff;
	if (reader.getElectronConductionCoeff(materialId, data))
	    electronConductionCoeff = BilinearInterpTable (spGrid, data);
	
	BilinearInterpTable ionConductionCoeff;
	if (reader.getIonConductionCoeff(materialId, data))
	    ionConductionCoeff = BilinearInterpTable (spGrid, data);
	
	BilinearInterpTable electronSpecificHeat;
	if (reader.getElectronSpecificHeat(materialId, data))
	    electronSpecificHeat = BilinearInterpTable(spGrid, data);
	
	BilinearInterpTable ionSpecificHeat;
	if (reader.getIonSpecificHeat(materialId, data))
	    ionSpecificHeat = BilinearInterpTable(spGrid, data);

	typedef MatTabMap::value_type MapPair;

	std::pair<MatTabMap::iterator, bool> retval;

	retval =
	    materials.insert(MapPair(materialId,
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

//---------------------------------------------------------------------------//
// has MaterialTables:
//   Return whether the material tables specified by the material id exists.
//---------------------------------------------------------------------------//

bool IMP::hasMaterialTables(int matId) const
{
    MatTabMap::const_iterator matit = materials.find(matId);

    return matit != materials.end();
}

//---------------------------------------------------------------------------//
// getMaterialTables:
//   Return the material tables specified by the material id.
//---------------------------------------------------------------------------//

const IMP::MaterialTables &IMP::getMaterialTables(int matId) const
{
    MatTabMap::const_iterator matit = materials.find(matId);

    // Make sure the material exists.
	
    Assert(matit != materials.end());

    return (*matit).second;
}

//---------------------------------------------------------------------------//
// getMaterialTables:
//   Return the material tables specified by the material id.
//---------------------------------------------------------------------------//

IMP::MaterialTables &IMP::getMaterialTables(int matId)
{
    MatTabMap::iterator matit = materials.find(matId);

    // Make sure the material exists.
	
    Assert(matit != materials.end());

    return (*matit).second;
}

//------------------------------------------------------------------------//
// interpolate:
//   Given a material state, a group number, and a pointer to a method
//   that returns a GroupedTable, return a unary operation on the
//   interpolated results from the table indicated by the method pointer
//   and group.
//------------------------------------------------------------------------//

template<class FT, class UnaryOperation>
void IMP::interpolate(const MaterialStateField<FT> &matState, int group,
		      PGroupedTable pTable, UnaryOperation op,
		      FT &results) const
{
    Require(&matState.matprops == this);
    Require(matState.size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < matState.size(); i++)
    {
	// Get the material tables for this location's particular material.
	
	const MaterialTables &mattabs = getMaterialTables(matState.getMatId(i));

	// Get the set of specific material property tables
	// (e.g. ionConduction tables).
	// (The set is over all of the groups.)

	const GroupedTable &groupedTable = (mattabs.*pTable)();
	
	// groups are 1-based

	Require(group >= 1 && group <= groupedTable.numGroups());

	// Get the specific table for this group number.
	
	const BilinearInterpTable &table = groupedTable.getTable(group);

	// Do the interpolation.

	typename FT::value_type intVal;

	if (table.hasData())
	    intVal = table.interpolate(matState.getMemento(i));
	else
	    intVal = 0.0;

	// Perform a function on the interpolated value.

	*resit++ = op(intVal);
    }
}

//---------------------------------------------------------------------------//
// interpolate:
//   Given a material state, and a pointer to a method
//   that returns a BilinearInterpTable, return a unary operation on the
//   interpolated results from the table indicated by the method pointer.
//---------------------------------------------------------------------------//

template<class FT, class UnaryOperation>
void IMP::interpolate(const MaterialStateField<FT> &matState,
		      PBilinearInterpTable pTable, UnaryOperation op,
		      FT &results) const
{
    Require(&matState.matprops == this);
    Require(matState.size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < matState.size(); i++)
    {
	// Get the material tables for this location's particular material.
	
	const MaterialTables &mattabs
	    = getMaterialTables(matState.getMatId(i));

	// Get the specific material property table (e.g. ionConduction table)
	
	const BilinearInterpTable &table = (mattabs.*pTable)();
	
	// Do the interpolation.

	typename FT::value_type intVal;

	if (table.hasData())
	    intVal = table.interpolate(matState.getMemento(i));
	else
	    intVal = 0.0;

	// Perform a function on the interpolated value.

	*resit++ = op(intVal);
    }
}

//---------------------------------------------------------------------------//
// getMaterialState:
//    Return a material state field from density, temperature, and materia id
//    fields.
//---------------------------------------------------------------------------//

template<class FT, class FT2>
IMP::MaterialStateField<FT> IMP::getMaterialState(const FT &density_,
						  const FT &electronTemp_,
						  const FT &ionTemp_,
						  const FT2 &matId_) const
{
    return MaterialStateField<FT>(*this, density_, electronTemp_,
				  ionTemp_, matId_);
}

//===========================================================================//
// InterpedMaterialProps::MaterialStateField<FT> Methods
//===========================================================================//

//---------------------------------------------------------------------------//
// InterpedMaterialProps::MaterialStateField<FT> Constructor:
//    Create an InterpedMaterialProps::MaterialStateField<FT>
//    density, temperature, and materia id fields.
//---------------------------------------------------------------------------//

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
	if (!matprops.hasMaterialTables(getMatId(i)))
	{
	    std::ostrstream os;
	    os << "InterpedMaterialProps::getMaterialState: "
	       << "Unable to find material "
	       << getMatId(i) << " in material table"
	       << " at " << i << "." << std::ends;
	    throw std::runtime_error(os.str());
	}
		
	const MaterialTables &mattabs = matprops.getMaterialTables(getMatId(i));
	
	const BilinearInterpGrid &grid = mattabs.getGrid();

	memento.push_back(grid.getMemento(getElectronTemp(i), getDensity(i)));
    }
}

//---------------------------------------------------------------------------//
// getElectronTemperature:
//    Return the temperature from the material state field.
//---------------------------------------------------------------------------//

template<class FT>
void IMP::MaterialStateField<FT>::getElectronTemperature(FT &results) const
{
    Require(size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < size(); i++)
	*resit++ = getElectronTemp(i);
}

//---------------------------------------------------------------------------//
// getIonTemperature:
//    Return the temperature from the material state field.
//---------------------------------------------------------------------------//

template<class FT>
void IMP::MaterialStateField<FT>::getIonTemperature(FT &results) const
{
    Require(size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < size(); i++)
	*resit++ = getIonTemp(i);
}

//---------------------------------------------------------------------------//
// getDensity:
//    Return the density from the material state field.
//---------------------------------------------------------------------------//

template<class FT>
void IMP::MaterialStateField<FT>::getDensity(FT &results) const
{
    Require(size() == results.size());
	
    FT::iterator resit = results.begin();
    
    for (int i=0; i < size(); i++)
	*resit++ = getDensity(i);
}

END_NS_XTM  // namespace XTM

//---------------------------------------------------------------------------//
//                              end of InterpedMaterialProps.cc
//---------------------------------------------------------------------------//
