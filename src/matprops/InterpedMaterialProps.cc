//----------------------------------*-C++-*----------------------------------//
// InterpedMaterialProps.cc
// Randy M. Roberts
// Wed Apr 15 08:44:49 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "InterpedMaterialProps.hh"

#include "BilinearInterpGrid.hh"
#include "BilinearInterpTable.hh"
#include "MaterialPropsReader.hh"

#include "ds++/SP.hh"
#include "ds++/Assert.hh"

#include <string>
using std::string;
#include <vector>
using std::vector;

#include <stdexcept>
#include <strstream>

typedef rtt_matprops::InterpedMaterialProps IMP;

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

	rtt_dsxx::SP<BilinearInterpGrid>
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

	rtt_dsxx::Mat2<double> data(tempGrid.size(), densityGrid.size(), 0.0);
	vector<BilinearInterpTable> tables(numGroups);

	// Read and load SigmaTotal interpolation table.
	
	for (int i=1; i<=numGroups; i++)
	{
	    if (reader.getSigmaTotal(materialId, i, data))
		tables[i-1] = BilinearInterpTable(spGrid, data);
	}

	GroupedTable sigmaTotal(energyUpperbounds, energyLowerbounds,
				tables);
	
	// Read and load SigmaAbsorption interpolation table.
	
	for (int i=1; i<=numGroups; i++)
	{
	    if (reader.getSigmaAbsorption(materialId, i, data))
		tables[i-1] = BilinearInterpTable(spGrid, data);
	}

	GroupedTable sigmaAbsorption(energyUpperbounds, energyLowerbounds,
				tables);
	
	// Read and load SigmaScattering interpolation table.
	
	for (int i=1; i<=numGroups; i++)
	{
	    if (reader.getSigmaScattering(materialId, i, data))
		tables[i-1] = BilinearInterpTable(spGrid, data);
	}

	GroupedTable sigmaScattering(energyUpperbounds, energyLowerbounds,
				tables);
	
	// Read and load SigmaEmission interpolation table.
	
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
						    sigmaScattering,
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

//---------------------------------------------------------------------------//
//                              end of InterpedMaterialProps.cc
//---------------------------------------------------------------------------//
