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
    Units units = Units::getAstroPhysUnits();

    std::ifstream ifs("testMatProp.inp");

    typedef FifiMatPropsReader::MaterialDefinition MatDef;
    vector<MatDef> matdefs;
    matdefs.push_back(MatDef("BeO, 13 Group", 2, 25.0));
    matdefs.push_back(MatDef("BeO, 3 Group", 3, 25.0));
    
    FifiMatPropsReader reader(matdefs, units, ifs);

    vector<int> matIds(matdefs.size());
    for (int i=0; i<matdefs.size(); i++)
	matIds[i] = matdefs[i].matid;
    
    IMP matProp(matIds, reader);

    cerr << "Made it past the matProp constructor" << endl;
    
    typedef vector<double> ccsf;

    const int ncells = 1;
    
    ccsf density(ncells);
    ccsf temp(ncells);
    vector<int> matid(ncells);

#if 0
    for (int i=0; i<ncells; i++)
    {
	density[i] = 0.2*i;
	temp[i] = 1.e-3*i;
	matid[i] = (i<5 ? matIds[0] : matIds[1]);
    }
#endif

    density[0] = 0.565;
    temp[0] = 1.0e+33;
    matid[0] = 3;
    
    cerr << "before getMaterialState" << endl;
    
    IMP::MaterialStateField<ccsf> matstate =
	matProp.getMaterialState(density, temp, temp, matid);

    cerr << "after getMaterialState" << endl;

    cout << std::scientific;

    cout << "Material State: " << endl;
    matstate.print(cout) << endl;
    
    ccsf results(ncells);
    
    matstate.getDensity(results);

    cout << "Density: " << endl;
    print(cout, results, 10);
    
    matstate.getElectronTemperature(results);

    cout << "Electron Temperature: " << endl;
    print(cout, results, 10);
    
    matstate.getSigmaTotal(1, results);

    cout << "SigmaTotal: " << endl;
    print(cout, results, 10);
    
    matstate.getSigmaAbsorption(1, results);

    cout << "SigmaAbsorption: " << endl;
    print(cout, results, 10);
    
    matstate.getElectronIonCoupling(results);

    cout << "ElectronIonCoupling: " << endl;
    print(cout, results, 10);
    
}
