//-------------------------------*-C++-*----------------------------------//
// FifiMatPropsReader.cc
// Randy M. Roberts
// Mon Apr 20 13:44:11 1998
//------------------------------------------------------------------------//
// @> 
//------------------------------------------------------------------------//

#include "FifiMatPropsReader.hh"
#include "radphys/RadiationPhysics.hh"
#include "units/Units.hh"

#include "ds++/Assert.hh"

#include "ds++/Mat.hh"
using rtt_dsxx::Mat2;

#include <iostream>
using std::istream;

#ifdef DEBUG_FIFIMATPROPSREADER
using std::cerr;
using std::endl;
#endif

#include <fstream>
using std::ifstream;

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <cctype>

#include <map>

#include <algorithm>

#include <strstream>
using std::istrstream;
using std::ostrstream;
using std::ends;

#include <cstdlib>

#include <stdexcept>

using namespace XTM;
using namespace rtt_matprops;

//------------------------------------------------------------------------//
// FMPR:
//    A shorthand typedef for the very long name, FifiMatPropsReader.
//------------------------------------------------------------------------//

typedef FifiMatPropsReader FMPR;

//------------------------------------------------------------------------//
// FifiMatPropsReader:
//   This constructor takes a vector of MaterialDefinition's,
//   the desired output Units, and the Fifi file's input stream.
//   The vector of MaterialDefinition's decides which materials to "look at"
//   from the file.  Any materials with ids not matching the "matid" field
//   of any of matdefs' entries will not be processed, and therefore, ignored.
//   The "abar" field of the MaterialDefinition's are the mass in amus of
//   the component material.
//------------------------------------------------------------------------//

FMPR::FifiMatPropsReader(const vector<MaterialDefinition> &matdefs,
			 const XTM::Units &outputUnits_, std::istream &is_)
    : MaterialPropsReader(outputUnits_), fifiParser(is_),
      fileUnits(XTM::Units::getAstroPhysUnits())
{
    // For fifi most of the units for the data are related
    // to Atronomical Physical Units.
    // Some will have to be converted to Astronomical Physical Units
    // before conversion to the user's output units.
    //
    // The file2OutputUnits object is the units converter from
    // file, i.e. Astronomical, units to the user's output units.
	
    file2OutputUnits = fileUnits / getOutputUnits();


    // Go throuth the supplied list of desired materials, and
    // process their data.
    
    for (vector<MaterialDefinition>::const_iterator mit = matdefs.begin();
	 mit != matdefs.end(); mit++)
    {
	const MaterialDefinition &mat = (*mit);
	int matid = mat.matid;

	// Make sure that there are no duplicates in the matdefs vector.
	
	if (hasMaterial(matid))
	{
	    ostrstream os;
	    os << "FifiMatPropsReader ctor: "
	       << "material: " << matid << " "
	       << " already processed." << ends;
	
	    throw std::runtime_error(os.str());
	}

	// If the Fifi file parser cannot find the material with this
	// matid, then this is a problem.
	
	if (!fifiParser.hasMaterial(matid))
	{
	    ostrstream os;
	    os << "FifiMatPropsReader ctor: "
	       << "material: " << matid << " "
	       << mat.name << " not found." << ends;
	
	    throw std::runtime_error(os.str());
	}

	// Create a nearly empty MaterialInfo.
	// It has no grids, only its matid and its mass in amus.
	// This MaterialInfo is then added to the map, where it
	// can be references by its matid, and subsequently filled with
	// its grid information.

	MaterialInfo matinfo(mat.name, matid, mat.abar);
	materialInfoMap.insert(std::make_pair(matid, matinfo));
    }

    // Go through all of the materials and fill in their temperature,
    // density, and group energy grid information.
    
    calcGridInfo();
}

//------------------------------------------------------------------------//
// getNextMaterial:
//   This method doesn't really do much but verify that the material exists,
//   and loads its name into the resulting string, and returns true.
//   If the material does not exist it returns false.
//------------------------------------------------------------------------//

bool FMPR::getNextMaterial(MaterialId materialId_, string &name)
{
    if (hasMaterial(materialId_))
    {
	const MaterialInfo &matinfo = getMaterialInfo(materialId_);
    	name = matinfo.name;
	return true;
    }
    return false;
}


//------------------------------------------------------------------------//
// calcGridInfo:
//    Go through all of the materials and fill in their temperature,
//    density, and group energy grid information.
//------------------------------------------------------------------------//

void FMPR::calcGridInfo()
{
    for (MatInfoMap::iterator matiter = materialInfoMap.begin();
	 matiter != materialInfoMap.end();
	 matiter++)
    {
	MaterialId matid = (*matiter).first;
	MaterialInfo &matInfo = (*matiter).second;

	// Get the temperature grid for this material

	fifiParser.getData(matid, "tgrid", matInfo.temperatureGrid);

	for (int i=0; i<matInfo.temperatureGrid.size(); i++)
	{
	    const XTM::Units &f2O = file2OutputUnits;

	    // The file units are in keV
	    // Convert these to user temperature units.
	    
	    matInfo.temperatureGrid[i] = f2O.ConvertTemperature(
		matInfo.temperatureGrid[i]);
	}
	
	// Get the density grid for this material
	
	fifiParser.getData(matid, "rgrid", matInfo.densityGrid);

	for (int i=0; i<matInfo.densityGrid.size(); i++)
	{
	    const XTM::Units &f2O = file2OutputUnits;

	    // The file units are in gm/cm**3
	    // Convert these to user density units.
	    
	    matInfo.densityGrid[i] = f2O.ConvertDensity(matInfo.densityGrid[i]);
	}
	
	// Get the energy grid for this material
	
	fifiParser.getData(matid, "hnugrid", matInfo.energyGrid);

	for (int i=0; i<matInfo.energyGrid.size(); i++)
	{
	    const XTM::Units &f2O = file2OutputUnits;

	    // The file units are in keV
	    // Convert these to user ***temperature*** units.
	    
	    matInfo.energyGrid[i] = f2O.ConvertTemperature(
		matInfo.energyGrid[i]);
	}

#ifdef DEBUG_FIFIMATPROPSREADER
	cerr << endl;
	cerr << "material: " << matid << " statistics:"
	     << " numTemperatures: " << matInfo.temperatureGrid.size()
	     << " numDensities: " << matInfo.densityGrid.size()
	     << " numGroups: " << matInfo.energyGrid.size() - 1
	     << endl;
#endif
    }
}

//------------------------------------------------------------------------//
// getMaterialInfo:
//   Return the MaterialInfo spec from the given materialId.
//   Throw an exception if the material is not in the map.
//------------------------------------------------------------------------//

const FMPR::MaterialInfo &FMPR::getMaterialInfo(MaterialId materialId) const
{
    MatInfoMap::const_iterator matiter = materialInfoMap.find(materialId);

    if (matiter == materialInfoMap.end())
    {
	ostrstream os;
	os << "FifiMatPropsReader::getMaterialInfo: "
	   << "material: " << materialId << " not found." << ends;
	
	throw std::runtime_error(os.str());
    }

    return (*matiter).second;    
}

//------------------------------------------------------------------------//
// getTemperatureGrid:
//   Return the temperature grid for the material given by materialId.
//   The call to getMaterialInfo will throw an exception if the material
//   is not found.
//------------------------------------------------------------------------//

bool FMPR::getTemperatureGrid(MaterialId materialId,
			      vector<double> &tempGrid_)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    tempGrid_ = matInfo.getTemperatureGrid();

    return true;
}

//------------------------------------------------------------------------//
// getDensityGrid:
//   Return the density grid for the material given by materialId.
//   The call to getMaterialInfo will throw an exception if the material
//   is not found.
//------------------------------------------------------------------------//

bool FMPR::getDensityGrid(MaterialId materialId,
			  vector<double> &densityGrid_)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    densityGrid_ = matInfo.getDensityGrid();

    return true;
}

//------------------------------------------------------------------------//
// getNumGroups:
//   Return the number of energy groups for the material given by materialId.
//   The call to getMaterialInfo will throw an exception if the material
//   is not found.
//------------------------------------------------------------------------//

bool FMPR::getNumGroups(MaterialId materialId, int &numGroups)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);
    numGroups = matInfo.getNumGroups();

    return true;
}

//------------------------------------------------------------------------//
// getEnergyUpperbounds:
//   Return the upper bounds of the specified group for the material
//   given by materialId.
//   The call to getMaterialInfo will throw an exception if the material
//   is not found.
//   Units: (temp)
//------------------------------------------------------------------------//

bool FMPR::getEnergyUpperbounds(MaterialId materialId, int group,
				double &energyUpperbounds_)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(group >= 1 && group <= matInfo.getNumGroups());

    // The upper bound for group ig is found at egrid[ig].
    // The lower bound for group ig is found at egrid[ig-1];
    
    energyUpperbounds_ = matInfo.getEnergyGrid()[group];

    return true;
}

//------------------------------------------------------------------------//
// getEnergyLowerbounds:
//   Return the lower bounds of the specified group for the material
//   given by materialId.
//   The call to getMaterialInfo will throw an exception if the material
//   is not found.
//   Units: (temp)
//------------------------------------------------------------------------//

bool FMPR::getEnergyLowerbounds(MaterialId materialId, int group,
				double &energyLowerbounds_)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(group >= 1 && group <= matInfo.getNumGroups());
    
    // The upper bound for group ig is found at egrid[ig].
    // The lower bound for group ig is found at egrid[ig-1];
    
    energyLowerbounds_ = matInfo.getEnergyGrid()[group-1];

    return true;
}

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

bool FMPR::getSigmaAbsorption(MaterialId materialId, int group,
			      Mat2<double> &dataMat)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(dataMat.nx() == matInfo.getNumTemperatures());
    Require(dataMat.ny() == matInfo.getNumDensities());

    if (fifiParser.hasKeyword(matInfo.matid, "ramg"))
    {
	getSigma(matInfo, group, "ramg", dataMat);
	return true;
    }

    return false;
}

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

bool FMPR::getSigmaScattering(MaterialId materialId, int group,
			      Mat2<double> &dataMat)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(dataMat.nx() == matInfo.getNumTemperatures());
    Require(dataMat.ny() == matInfo.getNumDensities());

    if (fifiParser.hasKeyword(matInfo.matid, "rsmg0"))
    {
	getSigma(matInfo, group, "rsmg0", dataMat);
	return true;
    }
    else if (fifiParser.hasKeyword(matInfo.matid, "rsmg"))
    {
	getSigma(matInfo, group, "rsmg", dataMat);
	return true;
    }

    return false;
}

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

bool FMPR::getSigmaTotal(MaterialId materialId, int group,
			 Mat2<double> &dataMat)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    int numTemps = matInfo.getNumTemperatures();
    int numDensities = matInfo.getNumDensities();

    Require(dataMat.nx() == matInfo.getNumTemperatures());
    Require(dataMat.ny() == matInfo.getNumDensities());

    // Does the parser have the absorption x-section and also
    // the l=0 scattering x-section?  If not just return false.
    
    if (!fifiParser.hasKeyword(materialId, "ramg") &&	
	!(fifiParser.hasKeyword(materialId, "rsmg")  ||
	  fifiParser.hasKeyword(materialId, "rsmg0")))
    {
	return false;
    }

    // Get the absorption x-section.
    // Have it initialized to zero, in case the fifi file does
    // not have an absorption x-section.
    
    Mat2<double> sigmaAbs(numTemps, numDensities, 0.0);

    if (fifiParser.hasKeyword(materialId, "ramg"))
	getSigma(matInfo, group, "ramg", sigmaAbs);
 
    // Get the l=0 scattering x-section.
    // Have it initialized to zero, in case the fifi file does
    // not have an l=0 scattering x-section.
    
    Mat2<double> sigmaSct(numTemps, numDensities, 0.0);

    if (fifiParser.hasKeyword(materialId, "rsmg0"))
    {
	getSigma(matInfo, group, "rsmg0", sigmaSct);
    }    
    else if (fifiParser.hasKeyword(materialId, "rsmg"))
    {
	getSigma(matInfo, group, "rsmg", sigmaSct);
    }

    // The total x-section is the sum of the absorption, and l=0
    // scattering x-sections.
    
    for (int i=0; i<numTemps; i++)
	for (int j=0; j<numDensities; j++)
	    dataMat(i,j) = sigmaAbs(i,j) + sigmaSct(i,j);

    return true;
}

//------------------------------------------------------------------------//
// getSigma:
//   Get a set of cross-section data, at the specified group for the
//   material given by the MaterialInfo spec, from the Fifi file parser
//   using the supplied keyword (e.g. "ramg", or "rsmg0").
//   The data are returned in a Mat2<double> dataMat(nTemps, nDens).
//   Units: (length^2/mass)
//------------------------------------------------------------------------//

void FMPR::getSigma(const MaterialInfo &matInfo, int group,
		    const string &keyword, Mat2<double> &dataMat)
{
    Require(dataMat.nx() == matInfo.getNumTemperatures());
    Require(dataMat.ny() == matInfo.getNumDensities());
    
    int numTemps = matInfo.getNumTemperatures();
    int numDensities = matInfo.getNumDensities();
    int numGroups = matInfo.getNumGroups();

    // The Fifi file parser will allocate the dataVec to what it
    // thinks is the correct size.  If we are really getting
    // cross-section data the size will be numGroups * numTemps * numDensities.
    
    vector<double> dataVec;
    fifiParser.getData(matInfo.matid, keyword, dataVec);

    Insist(dataVec.size() == numGroups * numTemps * numDensities,
	   (string("FifiMatPropsReader::getSigma: ") +
	    "data for keyword: \"" + keyword +
	    "\" is not number of groups by number of temperatures" +
	    " by number of densities.").c_str());
    
    vector<double>::const_iterator dataiter = dataVec.begin();

    // Strip out the data for the desired group.
    // Place the data into an numTemps by numDensities Mat2<double>.
    // Convert the x-sections to the desired output units.
    
    for (int i=0; i<numTemps; i++)
    {
	for (int j=0; j<numDensities; j++)
	{
	    const XTM::Units &f2O = file2OutputUnits;

	    const double datum = *(dataiter + group - 1);

	    // The file data is in units of cm**2/gm
	    // Convert this to the output units.
	    
	    dataMat(i,j) = f2O.InvertMass(f2O.ConvertLength(datum, 2));

	    dataiter += numGroups;
	}
    }
}

//------------------------------------------------------------------------//
// hasMaterial:
//    Returns true if material is found based on the mapping from
//    the materialId.
//------------------------------------------------------------------------//

bool FMPR::hasMaterial(MaterialId materialId) const
{
    return materialInfoMap.find(materialId) != materialInfoMap.end();
}

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

bool FMPR::getSigmaEmission(MaterialId materialId, int group,
			    Mat2<double> &data)
{
    // ASSUMES LTE!!!
    
    return getSigmaAbsorption(materialId, group, data);
}

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

bool FMPR::getElectronIonCoupling(MaterialId materialId, Mat2<double> &data)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(data.nx() == matInfo.getNumTemperatures());
    Require(data.ny() == matInfo.getNumDensities());
    
    const double abar = matInfo.abar;
    
    vector<double> z;
    if (!fifiParser.getData(matInfo.matid, "tfree", z))
	return false;

    // No units conversion necessary for "tfree" electrons/ion.

    // Create RadiationPhysics for the output units.

    RadiationPhysics radPhys(getOutputUnits());
    
    for (int iz=0, i=0; i<matInfo.getNumTemperatures(); i++)
	for (int j=0; j<matInfo.getNumDensities(); j++)
	    radPhys.getElectIonCoupling(matInfo.getDensityGrid()[j],
					matInfo.getTemperatureGrid()[i],
					z[iz++], abar, data(i,j));
    return true;
}

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

bool FMPR::getElectronConductionCoeff(MaterialId materialId, Mat2<double> &data)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(data.nx() == matInfo.getNumTemperatures());
    Require(data.ny() == matInfo.getNumDensities());
    
    const double abar = matInfo.abar;
    
    vector<double> z;
    if (!fifiParser.getData(matInfo.matid, "tfree", z))
	return false;

    // No units converselectron necessary for "tfree" electrons/ion.

    // Create RadiationPhysics for the output units.

    RadiationPhysics radPhys(getOutputUnits());
    
    for (int iz=0, i=0; i<matInfo.getNumTemperatures(); i++)
	for (int j=0; j<matInfo.getNumDensities(); j++)
	    radPhys.getElectronConductionCoeff(matInfo.getDensityGrid()[j],
					       matInfo.getTemperatureGrid()[i],
					       z[iz++], abar, data(i,j));
    return true;
}

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

bool FMPR::getIonConductionCoeff(MaterialId materialId, Mat2<double> &data)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(data.nx() == matInfo.getNumTemperatures());
    Require(data.ny() == matInfo.getNumDensities());
    
    const double abar = matInfo.abar;
    
    vector<double> z;
    if (!fifiParser.getData(matInfo.matid, "tfree", z))
	return false;

    // No units conversion necessary for "tfree" electrons/ion.

    // Create RadiationPhysics for the output units.

    RadiationPhysics radPhys(getOutputUnits());
    
    for (int iz=0, i=0; i<matInfo.getNumTemperatures(); i++)
	for (int j=0; j<matInfo.getNumDensities(); j++)
	    radPhys.getIonConductionCoeff(matInfo.getDensityGrid()[j],
					  matInfo.getTemperatureGrid()[i],
					  z[iz++], abar, data(i,j));
    return true;
}

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

bool FMPR::getElectronSpecificHeat(MaterialId materialId, Mat2<double> &data)
{
    // The electron specific heat is calculated as the temperature
    // derivative of the internal electron energy per unit mass.
    
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(data.nx() == matInfo.getNumTemperatures());
    Require(data.ny() == matInfo.getNumDensities());
    
    Mat2<double> eelect(matInfo.getNumTemperatures(),
			matInfo.getNumDensities());

    // Begin a scoping block.
    {
	vector<double> eelect_v;
	if (!fifiParser.getData(matInfo.matid, "eelect", eelect_v))
	    return false;

	// File units for internal energy are: MJ/kg (MegaJoules / KiloGram)
	const Units &OU = getOutputUnits();

	for (int iz=0, i=0; i<matInfo.getNumTemperatures(); i++)
	    for (int j=0; j<matInfo.getNumDensities(); j++)
		eelect(i, j) = 1.0e6 *
		    OU.InvertEnergy(OU.ConvertMass(eelect_v[iz++]));
    }

    calcTemperatureDerivative(materialId, eelect, data);

    return true;
}

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

bool FMPR::getIonSpecificHeat(MaterialId materialId, Mat2<double> &data)
{
    // The ion specific heat is calculated as the temperature
    // derivative of the internal ion energy per unit mass.
    
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(data.nx() == matInfo.getNumTemperatures());
    Require(data.ny() == matInfo.getNumDensities());
    
    Mat2<double> enuc(matInfo.getNumTemperatures(),
		      matInfo.getNumDensities());

    // Begin a scoping block.
    {
	vector<double> enuc_v;
	if (!fifiParser.getData(matInfo.matid, "enuc", enuc_v))
	    return false;

	// File units for internal energy are: MJ/kg (MegaJoules / KiloGram)
	
	const Units &OU = getOutputUnits();

	for (int iz=0, i=0; i<matInfo.getNumTemperatures(); i++)
	    for (int j=0; j<matInfo.getNumDensities(); j++)
		enuc(i, j) = 1.0e6 *
		    OU.InvertEnergy(OU.ConvertMass(enuc_v[iz++]));
    }
    
    calcTemperatureDerivative(materialId, enuc, data);

    return true;
}

//------------------------------------------------------------------------//
// calcTemperatureDerivative:
//    Calculate the temperature derivative of the nTemps x nDensities
//    quantities.
//------------------------------------------------------------------------//

void FMPR::calcTemperatureDerivative(MaterialId materialId,
				     const Mat2<double> &data,
				     Mat2<double> &derivative) const
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(data.nx() == matInfo.getNumTemperatures());
    Require(data.ny() == matInfo.getNumDensities());
    Require(derivative.nx() == matInfo.getNumTemperatures());
    Require(derivative.ny() == matInfo.getNumDensities());

    const vector<double> &tgrid = matInfo.temperatureGrid;

    int numT = matInfo.getNumTemperatures();
    
    for (int j=0; j<matInfo.getNumDensities(); j++)
    {
	// The first two points are used for the derivative's first point.
	
	derivative(0, j) = (data(1, j) - data(0,j)) /
	    (tgrid[1] - tgrid[0]);

	// We will approximate the function with a quadratic interpolation,
	// and take the first deriviative of the quadratic
	// as the derivative.
	
	for (int i=1; i<numT-1; i++)
	{
	    const double dy21 = data(i, j)  - data(i-1, j);
	    const double dy32 = data(i+1, j)- data(i, j);
	    const double dt32 = tgrid[i+1]-tgrid[i];
	    const double dt21 = tgrid[i]-tgrid[i-1];
	    const double dt31 = tgrid[i+1]-tgrid[i-1];

	    const double dydt = (dt32*dy21/dt21 + dt21*dy32/dt32) / dt31;
	    
	    derivative(i, j) = dydt;
	}

	// The last two points are used for the derivative's last point.
	
	derivative(numT-1, j) = (data(numT-1, j) - data(numT-2, j)) /
	    (tgrid[numT-1] - tgrid[numT-2]);
    }
}

//------------------------------------------------------------------------//
//                              end of FifiMatPropsReader.cc
//------------------------------------------------------------------------//
