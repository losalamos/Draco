#include "matprops/InterpedMaterialProps.hh"

#include "matprops/BilinearInterpGrid.hh"
#include "matprops/BilinearInterpTable.hh"
#include "matprops/test/TestMatPropReader.hh"

#include "3T/Units.hh"

#include "ds++/SP.hh"
#include "ds++/Assert.hh"

#include <string>
#include <vector>
using std::vector;
#include <iostream>
using std::cerr;
using std::endl;

using namespace XTM;

typedef InterpedMaterialProps IMP;

void testMatProp()
{
    Units units = Units::getSIUnits();
    
    TestMatPropReader reader(units);

    IMP matProp(reader);

    typedef vector<double> ccsf;

    const int ncells = 5;
    
    ccsf density(ncells);
    ccsf temp(ncells);
    vector<int> matid(ncells);

    for (int i=0; i<ncells; i++)
    {
	density.push_back(10.0*i);
	temp.push_back(i);
	matid.push_back(i<5 ? 1 : 10);
    }
    
    IMP::MaterialStateField<ccsf> matstate =
	matProp.getMaterialState(density, temp, temp, matid);

    ccsf results(ncells);
    
    matstate.getSigmaTotal(1, results);
    matstate.getElectronIonCoupling(results);
}
