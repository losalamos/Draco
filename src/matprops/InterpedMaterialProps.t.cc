//----------------------------------*-C++-*----------------------------------//
// InterpedMaterialProps.t.cc
// Randy M. Roberts
// Tue Sep 01 14:27:00 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/InterpedMaterialProps.hh"

#include "matprops/BilinearInterpGrid.hh"
#include "matprops/BilinearInterpTable.hh"
#include "matprops/MaterialPropsReader.hh"

#include "ds++/SP.hh"
#include "ds++/Assert.hh"

#include <stdexcept>
#include <strstream>

typedef rtt_matprops::InterpedMaterialProps IMP;

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
    Require(matState.pMatprops == this);
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
    Require(matState.pMatprops == this);
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
    : pMatprops(&matprops_), theSize(density_.size()),
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
	if (!pMatprops->hasMaterialTables(getMatId(i)))
	{
	    std::ostrstream os;
	    os << "InterpedMaterialProps::getMaterialState: "
	       << "Unable to find material "
	       << getMatId(i) << " in material table"
	       << " at " << i << "." << std::ends;
	    throw std::runtime_error(os.str());
	}
		
	const MaterialTables &mattabs =
	    pMatprops->getMaterialTables(getMatId(i));
	
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

