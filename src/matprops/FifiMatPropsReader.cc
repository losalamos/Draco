//----------------------------------*-C++-*----------------------------------//
// FifiMatPropsReader.cc
// Randy M. Roberts
// Mon Apr 20 13:44:11 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/FifiMatPropsReader.hh"
#include "radphys/RadiationPhysics.hh"
#include "units/Units.hh"

#include "ds++/Assert.hh"

#include "ds++/Mat.hh"
using dsxx::Mat2;

#include <iostream>
using std::istream;
using std::cin;
using std::cerr;
using std::cout;
using std::endl;

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

typedef FifiMatPropsReader FMPR;

FMPR::FifiMatPropsReader(const vector<MaterialDefinition> &matdefs,
			 const Units &outputUnits_, std::istream &is_)
    : MaterialPropsReader(outputUnits_), fifiParser(is_),
      fileUnits(Units::getAstroPhysUnits())
{
    // For fifi most of the units for the data are related
    // to Atronomical Physical Units.
    // Some will have to be converted to Astronomical Physical Units
    // before conversion to the user's output units.
    //
    // The file2OutputUnits object is the units converter from
    // file, i.e. Astronomical, units to the user's output units.
	
    file2OutputUnits = fileUnits / getOutputUnits();

    for (vector<MaterialDefinition>::const_iterator mit = matdefs.begin();
	 mit != matdefs.end(); mit++)
    {
	const MaterialDefinition &mat = (*mit);
	int matid = mat.matid;
	
	if (hasMaterial(matid))
	{
	    ostrstream os;
	    os << "FifiMatPropsReader ctor: "
	       << "material: " << matid << " "
	       << " already processed." << ends;
	
	    throw std::runtime_error(os.str());
	}
	
	if (!fifiParser.hasMaterial(matid))
	{
	    ostrstream os;
	    os << "FifiMatPropsReader ctor: "
	       << "material: " << matid << " "
	       << mat.name << " not found." << ends;
	
	    throw std::runtime_error(os.str());
	}

	MaterialInfo matinfo(matid, mat.abar);
	materialInfoMap.insert(std::make_pair(matid, matinfo));
    }
    
    calcGridInfo();
}

bool FMPR::getNextMaterial(MaterialId materialId_, string &name)
{
    if (hasMaterial(materialId_))
    {
	std::ostrstream os;
	os << "material" << materialId_ << std::ends;
	name = os.str();
	return true;
    }
    return false;
}


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

	cerr << "matInfo: " << matid << " tgrid: " << endl;
	for (int i=0; i<matInfo.temperatureGrid.size(); i++)
	{
	    const Units &f2O = file2OutputUnits;

	    // The file units are in keV
	    // Convert these to user temperature units.
	    
	    matInfo.temperatureGrid[i] = f2O.ConvertTemperature(
		matInfo.temperatureGrid[i]);
	    cerr << matInfo.temperatureGrid[i] << " ";
	}
	cerr << endl;
	
	// Get the density grid for this material
	
	fifiParser.getData(matid, "rgrid", matInfo.densityGrid);

	cerr << "matInfo: " << matid << " rgrid: " << endl;
	for (int i=0; i<matInfo.densityGrid.size(); i++)
	{
	    const Units &f2O = file2OutputUnits;

	    // The file units are in gm/cm**3
	    // Convert these to user density units.
	    
	    matInfo.densityGrid[i] = f2O.ConvertDensity(matInfo.densityGrid[i]);
	    cerr << matInfo.densityGrid[i] << " ";
	}
	cerr << endl;
	
	// Get the energy grid for this material
	
	fifiParser.getData(matid, "hnugrid", matInfo.energyGrid);

	for (int i=0; i<matInfo.energyGrid.size(); i++)
	{
	    const Units &f2O = file2OutputUnits;

	    // The file units are in keV
	    // Convert these to user ***temperature*** units.
	    
	    matInfo.energyGrid[i] = f2O.ConvertTemperature(
		matInfo.energyGrid[i]);
	}
	
	cerr << endl;
	cerr << "material: " << matid << " statistics:"
	     << " numTemperatures: " << matInfo.temperatureGrid.size()
	     << " numDensities: " << matInfo.densityGrid.size()
	     << " numGroups: " << matInfo.energyGrid.size() - 1
	     << endl;
    }
}

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

bool FMPR::getTemperatureGrid(MaterialId materialId,
			      vector<double> &tempGrid_)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    tempGrid_ = matInfo.getTemperatureGrid();

    return true;
}

bool FMPR::getDensityGrid(MaterialId materialId,
			  vector<double> &densityGrid_)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    densityGrid_ = matInfo.getDensityGrid();

    return true;
}

bool FMPR::getNumGroups(MaterialId materialId, int &numGroups)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);
    numGroups = matInfo.getNumGroups();

    return true;
}

bool FMPR::getEnergyUpperbounds(MaterialId materialId, int group,
				double &energyUpperbounds_)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(group >= 1 && group <= matInfo.getNumGroups());
    
    energyUpperbounds_ = matInfo.getEnergyGrid()[group];

    return true;
}

bool FMPR::getEnergyLowerbounds(MaterialId materialId, int group,
				double &energyLowerbounds_)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(group >= 1 && group <= matInfo.getNumGroups());
    
    energyLowerbounds_ = matInfo.getEnergyGrid()[group-1];

    return true;
}

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

bool FMPR::getSigmaTotal(MaterialId materialId, int group,
			 Mat2<double> &dataMat)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    int numTemps = matInfo.getNumTemperatures();
    int numDensities = matInfo.getNumDensities();

    Require(dataMat.nx() == matInfo.getNumTemperatures());
    Require(dataMat.ny() == matInfo.getNumDensities());
    
    if (!fifiParser.hasKeyword(materialId, "ramg") &&	
	!(fifiParser.hasKeyword(materialId, "rams")  ||
	  fifiParser.hasKeyword(materialId, "rams0")))
    {
	return false;
    }
	  
    Mat2<double> sigmaAbs(numTemps, numDensities, 0.0);

    if (fifiParser.hasKeyword(materialId, "ramg"))
	getSigma(matInfo, group, "ramg", sigmaAbs);
 
    Mat2<double> sigmaSct(numTemps, numDensities, 0.0);

    if (fifiParser.hasKeyword(materialId, "rsmg0"))
    {
	getSigma(matInfo, group, "rsmg0", sigmaSct);
    }    
    else if (fifiParser.hasKeyword(materialId, "rsmg"))
    {
	getSigma(matInfo, group, "rsmg", sigmaSct);
    }

    for (int i=0; i<numTemps; i++)
	for (int j=0; j<numDensities; j++)
	    dataMat(i,j) = sigmaAbs(i,j) + sigmaSct(i,j);

    return true;
}

void FMPR::getSigma(const MaterialInfo &matInfo, int group,
		    const string &keyword, Mat2<double> &dataMat)
{
    Require(dataMat.nx() == matInfo.getNumTemperatures());
    Require(dataMat.ny() == matInfo.getNumDensities());
    
    int numTemps = matInfo.getNumTemperatures();
    int numDensities = matInfo.getNumDensities();
    int numGroups = matInfo.getNumGroups();
    
    vector<double> dataVec;
    fifiParser.getData(matInfo.matid, keyword, dataVec);

    Insist(dataVec.size() == numGroups * numTemps * numDensities,
	   (string("FifiMatPropsReader::getSigma: ") +
	    "data for keyword: \"" + keyword +
	    "\" is not number of groups by number of temperatures" +
	    " by number of densities.").c_str());
    
    vector<double>::const_iterator dataiter = dataVec.begin();
    
    for (int i=0; i<numTemps; i++)
    {
	for (int j=0; j<numDensities; j++)
	{
	    const Units &f2O = file2OutputUnits;

	    const double datum = *(dataiter + group - 1);

	    // The file data is in units of cm**2/gm
	    // Convert this to the output units.
	    
	    dataMat(i,j) = f2O.InvertMass(f2O.ConvertLength(datum, 2));

	    dataiter += numGroups;
	}
    }
}

bool FMPR::hasMaterial(MaterialId materialId) const
{
    return materialInfoMap.find(materialId) != materialInfoMap.end();
}

bool FMPR::getSigmaEmission(MaterialId materialId, int group,
			    Mat2<double> &data)
{
    // ASSUMES LTE!!!
    
    return getSigmaAbsorption(materialId, group, data);
}

bool FMPR::getElectronIonCoupling(MaterialId materialId, Mat2<double> &data)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(data.nx() == matInfo.getNumTemperatures());
    Require(data.ny() == matInfo.getNumDensities());
    
    const double abar = matInfo.abar;
    
    vector<double> z;
    if (!fifiParser.getData(matInfo.matid, "tfree", z))
	return false;

    // No units conversion necessary for "tfree".

    // Create RadiationPhysics for the output units.

    RadiationPhysics radPhys(getOutputUnits());
    
    for (int iz=0, i=0; i<matInfo.getNumTemperatures(); i++)
	for (int j=0; j<matInfo.getNumDensities(); j++)
	    radPhys.getElectIonCoupling(matInfo.getDensityGrid()[j],
					matInfo.getTemperatureGrid()[i],
					z[iz++], abar, data(i,j));
    return true;
}

bool FMPR::getElectronConductionCoeff(MaterialId materialId, Mat2<double> &data)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(data.nx() == matInfo.getNumTemperatures());
    Require(data.ny() == matInfo.getNumDensities());
    
    const double abar = matInfo.abar;
    
    vector<double> z;
    if (!fifiParser.getData(matInfo.matid, "tfree", z))
	return false;

    // No units converselectron necessary for "tfree".

    // Create RadiationPhysics for the output units.

    RadiationPhysics radPhys(getOutputUnits());
    
    for (int iz=0, i=0; i<matInfo.getNumTemperatures(); i++)
	for (int j=0; j<matInfo.getNumDensities(); j++)
	    radPhys.getElectronConductionCoeff(matInfo.getDensityGrid()[j],
					       matInfo.getTemperatureGrid()[i],
					       z[iz++], abar, data(i,j));
    return true;
}

bool FMPR::getIonConductionCoeff(MaterialId materialId, Mat2<double> &data)
{
    const MaterialInfo &matInfo = getMaterialInfo(materialId);

    Require(data.nx() == matInfo.getNumTemperatures());
    Require(data.ny() == matInfo.getNumDensities());
    
    const double abar = matInfo.abar;
    
    vector<double> z;
    if (!fifiParser.getData(matInfo.matid, "tfree", z))
	return false;

    // No units conversion necessary for "tfree".

    // Create RadiationPhysics for the output units.

    RadiationPhysics radPhys(getOutputUnits());
    
    for (int iz=0, i=0; i<matInfo.getNumTemperatures(); i++)
	for (int j=0; j<matInfo.getNumDensities(); j++)
	    radPhys.getIonConductionCoeff(matInfo.getDensityGrid()[j],
					  matInfo.getTemperatureGrid()[i],
					  z[iz++], abar, data(i,j));
    return true;
}

bool FMPR::getElectronSpecificHeat(MaterialId materialId, Mat2<double> &data)
{
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

bool FMPR::getIonSpecificHeat(MaterialId materialId, Mat2<double> &data)
{
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

//---------------------------------------------------------------------------//
//                              end of FifiMatPropsReader.cc
//---------------------------------------------------------------------------//
