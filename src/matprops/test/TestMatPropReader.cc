//----------------------------------*-C++-*----------------------------------//
// TestMatPropReader.cc
// Randy M. Roberts
// Mon Apr 20 15:55:23 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "matprops/test/TestMatPropReader.hh"

#include "matprops/FifiMatPropsReader.hh"
#include "units/Units.hh"

#include <fstream>
#include <strstream>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <string>
using std::string;
#include <vector>
using std::vector;

using namespace rtt_matprops;

TestMatPropReader::TestMatPropReader(const XTM::Units &units,
				     const string &filename,
				     const vector<int> &materialIds)
{
    std::ifstream ifs(filename.c_str());
    FifiMatPropsReader reader(units, ifs);

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

	vector<double> data;

	reader.getTemperatureGrid(materialId, data);

	cout << "Temperature Grid for material: "
	     << materialId << endl;

	for (int j=0; j<data.size(); j++)
	    cout << " " << data[j];
	cout << endl;
    }
}

//---------------------------------------------------------------------------//
//                              end of TestMatPropReader.cc
//---------------------------------------------------------------------------//
