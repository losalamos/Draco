//----------------------------------*-C++-*----------------------------------//
// testMmcMatProp.cc
// John McGhee
// Mon Sep 21 09:53:46 1998
//---------------------------------------------------------------------------//
// @> Tests the Multi-Material Cell Material Properties Class.
//---------------------------------------------------------------------------//

#include <fstream>
#include <strstream>
#include "matprops/test/testMmcMatProp.hh"
#include <vector>
#include "matprops/FifiMatPropsReader.hh"
#include "ds++/SP.hh"
#include "matprops/InterpedMaterialProps.hh"
#include "matprops/MultiMatCellMatProps.hh"

using XTM::InterpedMaterialProps;
using XTM::FifiMatPropsReader;
using dsxx::SP;

testMmcMatProp::testMmcMatProp()
{
//empty
}

testMmcMatProp::~testMmcMatProp()
{
//empty
}

void testMmcMatProp::execute_test()
{

// Create a hypothetical set of problem parameters.


    int ncell  = 4;
    std::vector<std::vector<double > > volume_fraction(ncell);
    std::vector<std::vector<double > > density(ncell);
    std::vector<std::vector<double > > electron_temp(ncell);
    std::vector<std::vector<double > > ion_temp(ncell);
    std::vector<std::vector<int    > > mat_id(ncell);
    for (int icell = 0; icell<ncell; icell++)
    {
	int nmat = icell%2 + 2;
	volume_fraction[icell].resize(nmat);
	density[icell].resize(nmat);
	electron_temp[icell].resize(nmat);
	ion_temp[icell].resize(nmat);
	mat_id[icell].resize(nmat);
	for (int imat = 0; imat<nmat; imat++)
	{
	    volume_fraction[icell][imat] = (icell+1)*(imat+1);
	    density[icell][imat] = 2.*volume_fraction[icell][imat];
	    electron_temp[icell][imat] = 2.*volume_fraction[icell][imat];
	    ion_temp[icell][imat] = 3.*volume_fraction[icell][imat];
	    mat_id[icell][imat] = imat+1;
	}
    }

// Create Material Properties Objects

    int nmat_prob = 3;
    std::vector<int> matId(nmat_prob);


    std::vector<FifiMatPropsReader::MaterialDefinition> mats(nmat_prob);
    for (int imat = 0; imat<nmat_prob; imat++)
    { 
	std::ostrstream os;
	os << "Material" << imat+1 << std::ends;
	double abar = 10.0*(imat+1);
	mats[imat] = FifiMatPropsReader::MaterialDefinition(os.str(),
							    imat+1, abar);
	matId[imat] = imat+1;
    }
    std::ifstream fifis("testMatProp.inp");
    FifiMatPropsReader reader(mats, XTM::Units(), fifis);
    SP<InterpedMaterialProps> spumcmp(new InterpedMaterialProps(matId,reader));

    typedef MultiMatCellMatProps<InterpedMaterialProps> MMCMP;

    MMCMP mmcmp(spumcmp);

// Create the material state field

    typedef std::vector<double> FT;
    typedef std::vector<std::vector<double> > FT1;
    typedef std::vector<std::vector<int> > FT2;

    MMCMP::MaterialStateField<FT,FT1> msf =
	mmcmp.getMaterialState<FT,FT1,FT2>(density, electron_temp,
			       ion_temp, volume_fraction, 
			       mat_id);
    
    
// Get some properties for examination

    {
	std::vector<double> results(ncell);
	std::vector<double> answer(ncell);
	int group = 1;
	msf.getSigmaAbsorption(group, results);
	for( int icell = 0; icell<ncell; icell++)
	{
	    std::cout << "results[" << icell << "] = " << results[icell] << std::endl;
	}
    }
	

}   

//
//---------------------------------------------------------------------------//
//                              end of testMmcMatProp.cc
//---------------------------------------------------------------------------//
