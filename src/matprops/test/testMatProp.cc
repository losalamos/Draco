#include "matprops/InterpedMaterialProps.hh"

#include "matprops/BilinearInterpGrid.hh"
#include "matprops/BilinearInterpTable.hh"
#include "matprops/FifiMatPropsReader.hh"

#include "units/Units.hh"

#include "ds++/SP.hh"
#include "ds++/Assert.hh"

#include <string>
#include <vector>
using std::vector;
#include <fstream>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <iomanip>
using std::scientific;
using std::setprecision;

using namespace XTM;

template<class T>
inline void print(std::ostream &os, T &v, int npline)
{
    int i = 0;
    for (T::const_iterator it = v.begin(); it != v.end(); it++)
    {
	os << *it << " ";
	if (!(++i % npline))
	    os << std::endl;
    }
    if (i%npline)
	os << std::endl;
}

typedef InterpedMaterialProps IMP;

void testMatProp()
{
    typedef vector<double> ccsf;
    typedef vector<int> ccif;

    Units units = Units::getAstroPhysUnits();
    // Units units;

    std::ifstream ifs("testMatProp.inp");

    typedef FifiMatPropsReader::MaterialDefinition MatDef;
    vector<MatDef> matdefs;
    matdefs.push_back(MatDef("BeO, 3 Group", 1, 25.0));
    matdefs.push_back(MatDef("BeO, 13 Group", 2, 25.0));
    matdefs.push_back(MatDef("BeO, 3 Group", 3, 25.0));
    
    FifiMatPropsReader reader(matdefs, units, ifs);

    vector<int> matIds(matdefs.size());
    for (int i=0; i<matdefs.size(); i++)
	matIds[i] = matdefs[i].matid;
    
    IMP matProp(matIds, reader);

    //.... DUMP SOME RESULTS TO SCREEN

    for (int i=0; i<matdefs.size(); i++)
    {
	int matid = matdefs[i].matid;
	cout << " mat#   = " << matid << endl;
	cout << " name   = " << matProp.getMaterialName(matid) << endl;
	cout << " info   = " << endl;
	cout << " ngroups= " << matProp.getNumGroups(matid) << endl;
	cout << " ndens  = " << matProp.getNumDensities(matid) << endl;
	cout << " ntemps = " << matProp.getNumTemperatures(matid) << endl;
	cout << " nscxs  = " << matProp.getMaxScatteringPnOrder(matid) << endl;

	vector<double> grid;

	grid = matProp.getTemperatureGrid(matid);
	ccsf lastTemp(1, grid[grid.size()-1]);
	
	cout << " Temperature Grid:" << endl;

	for (int j = 0; j < grid.size(); )
	{
	    cout << " ";
	    for (int jj=0; jj < 6 && j < grid.size() ; jj++, j++)
		cout << " " << scientific << setprecision(4) << grid[j];
	    cout << endl;
	}

	grid = matProp.getDensityGrid(matid);
	ccsf lastDens(1,grid[grid.size()-1]);

	cout << " Density Grid:" << endl;
	
	for (int j = 0; j < grid.size(); )
	{
	    cout << " ";
	    for (int jj=0; jj < 6 && j < grid.size() ; jj++, j++)
		cout << " " << scientific << setprecision(4) << grid[j];
	    cout << endl;
	}

	ccif matids(1);
	matids[0] = matid;
	
	IMP::MaterialStateField<ccsf> matstate =
	    matProp.getMaterialState(lastDens, lastTemp, lastTemp, matids);
	
	cout << " ABS(ndens,ntemps,1:ngroups):" << endl;
	for (int j = 1; j <= matProp.getNumGroups(matid); )
	{
	    cout << " ";
	    for (int jj=0; jj < 6 && j <= matProp.getNumGroups(matid);
		 jj++, j++)
	    {
		ccsf abs(1);
		matstate.getSigmaAbsorption(j, abs);
		cout << " " << scientific << setprecision(4)
		     << abs[0] / lastDens[0];
	    }
	    cout << endl;
	}

	cout << endl;
    }
    
    const int ncells = 200;
    int matnum = 0;
    int ngroups = matProp.getNumGroups(matdefs[matnum].matid);
    
    ccsf density(ncells);
    ccsf temp(ncells);
    vector<int> matid(ncells);
    vector<ccsf> abs(ngroups, ccsf(ncells));
    vector<ccsf> sigt(ngroups, ccsf(ncells));
    vector<ccsf> emission(ngroups, ccsf(ncells));
    ccsf cve(ncells);
    ccsf cvi(ncells);
    ccsf eic(ncells);
    ccsf tce(ncells);
    ccsf tci(ncells);

#if 0
    for (int i=0; i<ngroups; i++)
    {
	abs[i] = ccsf(ncells);
	sigt[i] = ccsf(ncells);
	emission[i] = ccsf(ncells);
    }
#endif

    for (int i=0; i<ncells; i++)
    {
	density[i] = units.InvertDensity((i+1)*1.5);
	matid[i] = matdefs[matnum].matid;
    }

    temp[0] = units.InvertTemperature(1.e4);
    for (int i=1; i<ncells; i++)
	temp[i] = 1.055*temp[i-1];
    
    IMP::MaterialStateField<ccsf> matstate =
	matProp.getMaterialState(density, temp, temp, matid);

    for (int ig = 1; ig <= ngroups; ig++)
    {
	matstate.getSigmaAbsorption(ig, abs[ig-1]);
	matstate.getSigmaTotal(ig, sigt[ig-1]);
	matstate.getSigmaEmission(ig, emission[ig-1]);
    }
    matstate.getElectronSpecificHeat(cve);
    matstate.getIonSpecificHeat(cvi);
    matstate.getElectronIonCoupling(eic);
    matstate.getElectronConductionCoeff(tce);
    matstate.getIonConductionCoeff(tci);

    int ig = ngroups;
    int m = 19;

    cout << scientific << setprecision(12);
    
    cout << " Material Properties Summary" << endl;
    cout << " cell# (m)==>" << m <<  "   group# (ig)==>" << ig << endl;
    cout << " Material Table# ==>" <<  matid[m-1] << endl;
    cout << endl;
    cout << " temp(m)     = " <<  temp[m-1] << endl;
    cout << " dens(m)     = " <<  density[m-1] << endl;
    cout << " abs(m,ig)   = " <<  abs[ig-1][m-1] << endl;
    cout << " sct(m,0,ig) = " <<  "no-op" << endl;
    cout << " sct(m,1,ig) = " <<  "no-op" << endl;
    cout << " sigt(m,ig)  = " <<  sigt[ig-1][m-1] << endl;
    cout << " rmv(m,ig)   = " <<  "no-op" << endl;
    cout << " emis(m,ig)  = " <<  emission[ig-1][m-1] << endl;
    cout << " cve(m)      = " <<  cve[m-1] << endl;
    cout << " cvi(m)      = " <<  cvi[m-1] << endl;
    cout << " eic(m)      = " <<  eic[m-1] << endl;
    cout << " tce(m)      = " <<  tce[m-1] << endl;
    cout << " tci(m)      = " <<  tci[m-1] << endl;
    cout << endl;
}
