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
using std::cerr;
using std::endl;

using namespace XTM;

typedef InterpedMaterialProps IMP;

void testMatProp()
{
    Units units = Units::getSIUnits();

    const int matIdsC[] = {1,2,3};
    vector<int> matIds(matIdsC, matIdsC+sizeof(matIdsC)/sizeof(int));

    std::ifstream ifs("testMatProp.inp");
    
    FifiMatPropsReader reader(units, ifs);

    IMP matProp(matIds, reader);

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
