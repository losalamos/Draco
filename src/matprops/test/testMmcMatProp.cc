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
#include <list>
#include <vector>
#include "matprops/FifiMatPropsReader.hh"
#include "ds++/SP.hh"
#include "matprops/InterpedMaterialProps.hh"
#include "matprops/MultiMatCellMatProps.hh"
#include <limits>

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

// Create a hypothetical set of problem parameters. Use lists as
// containers to stress the system, as list has no [] operator.
// This data is used with the fifi opacity input file to define the
// problem.

    using std::cout;
    using std::endl;
    using std::list;

    int ncell  = 4;
    typedef list<double >        FTVD;
    typedef list<list<double > > FTVVD;
    typedef list<list<int    > > FTVVI;

    FTVVD volume_fraction;
    FTVVD density;
    FTVVD electron_temp;
    FTVVD ion_temp;
    FTVVI mat_id;

    for (int icell = 0; icell<ncell; icell++)
    {
	int nmat = icell%2 + 2;

	volume_fraction.push_back(list<double>());
	density.push_back(list<double>());
	electron_temp.push_back(list<double>());
	ion_temp.push_back(list<double>());
	mat_id.push_back(list<int>());

	for (int imat = 0; imat<nmat; imat++)
	{
	    volume_fraction.back().push_back((icell+1)*(imat+1));
	    density.back().push_back(2.125*volume_fraction.back().back());
	    electron_temp.back().push_back(3.52*volume_fraction.back().back());
	    ion_temp.back().push_back(3.771*volume_fraction.back().back());
	    mat_id.back().push_back(imat+1);
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
    std::ifstream fifis("testMmcMatProp.opac.inp");
    FifiMatPropsReader reader(mats, XTM::Units(), fifis);
    SP<InterpedMaterialProps> spumcmp(new InterpedMaterialProps(matId,reader));

    typedef MultiMatCellMatProps<InterpedMaterialProps> MMCMP;

    MMCMP mmcmp(spumcmp);

// Create the material state field

    typedef FTVD   FT;
    typedef FTVVD  FT1;
    typedef FTVVI  FT2;

    MMCMP::MaterialStateField<FT,FT1> msf =
	mmcmp.getMaterialState<FT,FT1,FT2>(density, electron_temp,
			       ion_temp, volume_fraction, 
			       mat_id);


	
// Print header.

    std::cout << std::endl;
    std::cout << std::scientific;
    std::cout.precision(5);
    std::cout << " --- Begin Multi-Material Cell Tests ---" << std::endl;    
    std::cout << std::endl;

// Define logical for self-test
    
    bool passed = true;
    double eps = 5000.*std::numeric_limits<double>::epsilon();
    std::cout << " Max Allowed Relative Error: " << eps 
	      << std::endl << std::endl;

// Echo Input

    std::cout.precision(5);
    std::cout << " ** Input Data **" << std::endl;
    std::cout << " Cell Mat Mat   Volume       Density  " << 
	"    Elec-Temp    Ion-Temp "
	      << std::endl;
    std::cout << "          ID#   (m**3)       (kg/m**3)" << 
	"       (K)         (K)"
	      << std::endl;
    
    FTVVD::iterator vit  = volume_fraction.begin();
    FTVVD::iterator dit  = density.begin();
    FTVVD::iterator teit = electron_temp.begin();
    FTVVD::iterator tiit = ion_temp.begin();
    FTVVI::iterator mit  = mat_id.begin();

    for( int icell = 0; icell<ncell; icell++)
    {
	int nmat = (*mit).size();

	FTVVD::value_type::iterator vitit  = (*vit++).begin();
	FTVVD::value_type::iterator ditit  = (*dit++).begin();
	FTVVD::value_type::iterator teitit = (*teit++).begin();
	FTVVD::value_type::iterator tiitit = (*tiit++).begin();
	FTVVI::value_type::iterator mitit  = (*mit++).begin();

	for (int imat = 0; imat<nmat; imat++)
	{
	    std::cout << " [" << icell << 
		"]  [" << imat << "]  " << 
		*mitit++  << "    " <<
		*vitit++  << "  " <<
		*ditit++  << "  " <<
		*teitit++ << "  " <<
		*tiitit++ << "  " << std::endl;
	}
	std::cout << std::endl;	
    }
    std::cout << std::endl;
    std::cout.precision(12);

// Check absorption opacity results.

    {
	std::cout << " ** Absorption Opacity **" << std::endl;
	std::cout << " Cell  Opacity (1/m)      Relative Error" 
		  << std::endl;
	std::list<double> results(ncell);
	std::vector<double> answer(ncell);
	std::vector<double> relErr(ncell);
	double relErrMax = 0.;
	answer[0] = 1.012916666667e+01;
	answer[1] = 5.475416666667e+01;
	answer[2] = 3.038750000000e+01;
	answer[3] = 1.095083333333e+02;
	int group = 1;
	msf.getSigmaAbsorption(group, results);
	FTVD::iterator resit = results.begin();
	for( int icell = 0; icell<ncell; icell++)
	{
	    relErr[icell] = std::abs(*resit-answer[icell])/answer[icell];
	    if (relErr[icell] > relErrMax) relErrMax = relErr[icell];
	    std::cout << " [" << icell << "]   " << 
		*resit++ << " " << relErr[icell] << std::endl;
	}
	std::cout << std::endl;
	passed = passed && relErrMax <= eps;
    }

// Check electron temperature results.

    {
	std::cout << " ** Electron Temperature **" << std::endl;
	std::cout << " Cell  Temperature (K)    Relative Error" 
		  << std::endl;
	std::list<double> results(ncell);
	std::vector<double> answer(ncell);
	std::vector<double> relErr(ncell);
	double relErrMax = 0.;
	answer[0] = 6.025766197348e+00;
	answer[1] = 1.932154402511e+01;
	answer[2] = 1.807729859204e+01;
    	answer[3] = 3.864308805021e+01;
	msf.getElectronTemperature(results);
	FTVD::iterator resit = results.begin();
	for( int icell = 0; icell<ncell; icell++)
	{
	    relErr[icell] = std::abs(*resit-answer[icell])/answer[icell];
	    if (relErr[icell] > relErrMax) relErrMax = relErr[icell];
	    std::cout << " [" << icell << "]   " << 
		*resit++ << " " << relErr[icell] << std::endl;
	}
	std::cout << std::endl;
	passed = passed && relErrMax <= eps;
    }
    
// Check electron specific heat results

    {
	std::cout << " ** Electron Specific Heat By Material**" << std::endl;
	std::cout << " Cell Mat   Cve (J/m**3-K)     Relative Error" 
		  << std::endl;
        std::list<std::list<double> > results(ncell);
        std::vector<std::vector<double> > answer(ncell);
	std::vector<std::vector<double> > relErr(ncell);
	double relErrMax = 0.; 
	FTVVD::iterator resit=results.begin();
	FTVVI::iterator mit= mat_id.begin();
	for (int icell = 0; icell<ncell; icell++)
	{
	    int nmat = (*mit++).size();
	    (*resit++).resize(nmat);
	    answer[icell].resize(nmat);
	    relErr[icell].resize(nmat);
	}
	answer[0] [0]= 7.225128263753e+00;
	answer[0] [1]= 8.925201530194e+00;
	answer[1] [0]= 1.445025652751e+01;
	answer[1] [1]= 1.785040306039e+01;
	answer[1] [2]= 6.757626427344e+01;
	answer[2] [0]= 2.167538479126e+01;
	answer[2] [1]= 2.677560459058e+01;
	answer[3] [0]= 2.890051305501e+01;
	answer[3] [1]= 3.570080612078e+01;
    	answer[3] [2]= 1.351525285469e+02;

	msf.getElectronSpecificHeatByMat(results);
	resit = results.begin();
	mit   = mat_id.begin();
	for( int icell = 0; icell<ncell; icell++)
	{
	    int nmat = (*mit++).size();
	    FTVD::iterator resitit= (*resit++).begin();
	    for (int imat = 0; imat<nmat; imat++)
	    {
		relErr[icell][imat] = 
		    std::abs(*resitit-answer[icell][imat])/
		    answer[icell][imat];
		if (relErr[icell][imat] > relErrMax)
		    relErrMax = relErr[icell][imat];
		std::cout << " [" << icell << 
		    "]  [" << imat << "]   " << 
		    *resitit++ << " " << relErr[icell][imat] << std::endl;
	    }
	    std::cout << std::endl;
	}
	std::cout << std::endl;
	passed = passed && relErrMax <= eps;
    }

// Print the status of the test.

    cout << endl;
    cout <<     " *********************************************" << endl;
    if (passed) 
    {
	cout << " ****   MMC Mat-Props Self Test: PASSED   ****" << endl;
    }
    else
    {
	cout << " ****   MMC Mat-Props Self Test: FAILED   ****" << endl;
    }
    cout <<     " *********************************************" << endl;
    cout << endl;

// Print footer.

    std::cout << std::endl;
    std::cout << " --- End Multi-Material Cell Tests ---" << std::endl;    
    std::cout << std::endl;

}   

//
//---------------------------------------------------------------------------//
//                              end of testMmcMatProp.cc
//---------------------------------------------------------------------------//
