//----------------------------------*-C++-*----------------------------------//
// testMmcMatProp.cc
// John McGhee
// Mon Sep 21 09:53:46 1998
//---------------------------------------------------------------------------//
// @> Tests the Multi-Material Cell Material Properties Class.
//---------------------------------------------------------------------------//

#include <sstream>
#include "testMmcMatProp.hh"
#include <list>
#include <vector>
#include "../FifiMatPropsReader.hh"
#include "ds++/SP.hh"
#include "../InterpedMaterialProps.hh"
#include "../MultiMatCellMatProps.hh"
#include <limits>

using rtt_matprops::MultiMatCellMatProps;
using rtt_matprops::InterpedMaterialProps;
using rtt_matprops::FifiMatPropsReader;
using rtt_dsxx::SP;

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
    std::vector<int> nmat(ncell);

    typedef list<double>     FTVD;
    typedef list<FTVD>       FTVVD;
    typedef list<list<int> > FTVVI;

    FTVVD volume_fraction;
    FTVVD density;
    FTVVD electron_temp;
    FTVVD ion_temp;
    FTVVI mat_id;
    
    
    for (int icell = 0; icell<ncell; icell++)
    {
	nmat[icell] = icell%2 + 2;

	volume_fraction.push_back(FTVD());
	density.push_back(FTVD());
	electron_temp.push_back(FTVD());
	ion_temp.push_back(FTVD());
	mat_id.push_back(list<int>());

	for (int imat = 0; imat<nmat[icell]; imat++)
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
	std::ostringstream os;
	os << "Material" << imat+1 << std::ends;
	double abar = 10.0*(imat+1);
	mats[imat] = FifiMatPropsReader::MaterialDefinition(os.str(),
							    imat+1, abar);
	matId[imat] = imat+1;
    }
    FifiMatPropsReader reader(mats, rtt_units::Units(),
                              "testMmcMatProp.opac.inp");
    SP<InterpedMaterialProps> spumcmp(new InterpedMaterialProps(matId,reader));

    typedef MultiMatCellMatProps<InterpedMaterialProps> MMCMP;

    MMCMP mmcmp(spumcmp);

// Create the material state field

    MMCMP::MaterialStateField<FTVD,FTVVD,FTVVI> msf =
	mmcmp.getMaterialState<FTVD,FTVVD,FTVVI>(density, electron_temp,
			       ion_temp, volume_fraction, 
			       mat_id);

// Print header.

    std::cout << std::scientific;
    std::cout << std::endl;
    std::cout << " --- Begin Multi-Material Cell Tests ---" << std::endl;    
    std::cout << std::endl;

// Define logical for self-test.  Define pass criteria.
    
    bool passed = true;
    bool results_OK;
    double eps = 5000.*std::numeric_limits<double>::epsilon();
    std::cout.precision(5);
    std::cout << " Max Allowed Relative Error: " << eps 
	      << std::endl << std::endl;

// Echo Input

    std::cout.precision(5);
    std::cout << " ** Cell Input Data **" << std::endl;
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

	FTVVD::value_type::iterator vitit  = (*vit++).begin();
	FTVVD::value_type::iterator ditit  = (*dit++).begin();
	FTVVD::value_type::iterator teitit = (*teit++).begin();
	FTVVD::value_type::iterator tiitit = (*tiit++).begin();
	FTVVI::value_type::iterator mitit  = (*mit++).begin();

	for (int imat = 0; imat<nmat[icell]; imat++)
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

    std::cout << " ** Fifi File Data ** (No temp or dens dependence):" 
	      << std::endl;
    std::cout << " Mat-ID  abs      Cve      Cvi      abar  tfree" << std::endl;
    std::cout << "       (m**2/kg) (J/kg-K) (J/kg-K)  amu " << std::endl;
    std::cout << "  1      2.3      3.4      3.6      10.0  3.0" << std::endl;
    std::cout << "  2      3.0      2.1      1.2      20.0  4.0" << std::endl;
    std::cout << "  3      7.0      5.3     10.1      30.0  5.0" << std::endl;
    std::cout << std::endl;

    std::cout.precision(12);

// Check the material temperature mapper

    {
	std::cout << " ** Electron Temperature Map **" << std::endl;
	std::cout << " Cell Mat   Temperature (K)    Relative Error" 
		  << std::endl;
	FTVD avgTemps(ncell);
	FTVD::iterator atit=avgTemps.begin();
	*atit++ = 6.025766197348e+00;
	*atit++ = 1.932154402511e+01;
	*atit++ = 1.807729859204e+01;
	*atit++ = 3.864308805021e+01;
        FTVVD results(ncell);
        std::vector<std::vector<double> > answer(ncell);
	FTVVD::iterator resit=results.begin();
	for (int icell = 0; icell<ncell; icell++)
	{
	    //	    (*resit++).resize(nmat[icell]);
	    answer[icell].resize(nmat[icell]);
	}
	answer[0]  [0]= 6.025766197348e+00;  
	answer[0]  [1]= 6.025766197348e+00;  
	answer[1]  [0]= 1.932154402511e+01;  
	answer[1]  [1]= 1.932154402511e+01; 
	answer[1]  [2]= 1.932154402511e+01; 
	answer[2]  [0]= 1.807729859204e+01; 
	answer[2]  [1]= 1.807729859204e+01; 
	answer[3]  [0]= 3.864308805021e+01;  
	answer[3]  [1]= 3.864308805021e+01; 
	answer[3]  [2]= 3.864308805021e+01; 
	msf.mapAvgElectronTemp(results, avgTemps);
	ByMatResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK;
    }



// Check electron energy deposition by material results.

    {
	std::cout << " ** Electron Energy Deposition By Material **" << 
	    std::endl;
	std::cout << " Cell Mat   Energy Dep (J)     Relative Error" 
		  << std::endl;
        FTVVD results(ncell);
	FTVD  newTemp(ncell);
	FTVD  volCell(ncell);
        std::vector<std::vector<double> > answer(ncell);
	FTVVD::iterator resit=results.begin();
	for (int icell = 0; icell<ncell; icell++)
	{
	    //	    (*resit++).resize(nmat[icell]);
	    answer[icell].resize(nmat[icell]);
	}
	answer[0]  [0]=   3.959370288537e+01;
	answer[0]  [1]=   3.498678999836e+01;
	answer[1]  [0]=   4.468019318305e+02;
	answer[1]  [1]=   6.012015750739e+02;
	answer[1]  [2]=   5.595314681841e+02;
	answer[2]  [0]=   1.589239212895e+03;
	answer[2]  [1]=   2.229872350304e+03;
	answer[3]  [0]=   1.936141704599e+03;
	answer[3]  [1]=   2.605206825320e+03;
	answer[3]  [2]=   2.424636362131e+03;
	FTVD::iterator ntit = newTemp.begin();
	*ntit++ = 9.0;
	*ntit++ = 22.5;
        *ntit++ = 35.0;
	*ntit++ = 45.0;
	FTVD::iterator volit = volCell.begin();
	*volit++ = 3.0;
	*volit++ = 12.0;
	*volit++ = 9.0;
	*volit++ = 13.0;
	msf.getElectronEnergyDepbyMat(results, newTemp, volCell);
	ByMatResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK;
    }

// Check ion energy deposition by material results.

    {
	std::cout << " ** Ion Energy Deposition By Material **" << 
	    std::endl;
	std::cout << " Cell Mat   Energy Dep (J)     Relative Error" 
		  << std::endl;
        FTVVD results(ncell);
	FTVD  newTemp(ncell);
	FTVD  volCell(ncell);
        std::vector<std::vector<double> > answer(ncell);
	FTVVD::iterator resit=results.begin();
	for (int icell = 0; icell<ncell; icell++)
	{
	    //	    (*resit++).resize(nmat[icell]);
	    answer[icell].resize(nmat[icell]);
	}
	answer[0]  [0]=   4.000261646872e+01;
	answer[0]  [1]=   1.487224093482e+01;
	answer[1]  [0]=   4.883241565177e+02;
	answer[1]  [1]=   3.433875986488e+02;
	answer[1]  [2]=   6.753223230259e+02;
	answer[2]  [0]=   1.630881198442e+03;
	answer[2]  [1]=   1.135982156343e+03;
	answer[3]  [0]=   2.049770074547e+03;
	answer[3]  [1]=   1.399609117623e+03;
	answer[3]  [2]=   1.252256725611e+03;
	FTVD::iterator ntit = newTemp.begin();
	*ntit++ = 9.0;
	*ntit++ = 23.5;
        *ntit++ = 35.0;
	*ntit++ = 46.0;
	FTVD::iterator volit = volCell.begin();
	*volit++ = 3.0;
	*volit++ = 12.0;
	*volit++ = 9.0;
	*volit++ = 13.0;
	msf.getIonEnergyDepbyMat(results, newTemp, volCell);
	ByMatResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK;
    }

// Check absorption opacity results.

    {
	std::cout << " ** Absorption Opacity **" << std::endl;
	std::cout << " Cell  Opacity (1/m)      Relative Error" 
		  << std::endl;
	std::vector<double> answer(ncell);
	answer[0] = 1.012916666667e+01;
	answer[1] = 5.475416666667e+01;
	answer[2] = 3.038750000000e+01;
	answer[3] = 1.095083333333e+02;
	FTVD results(ncell);
	int group = 1;
	msf.getSigmaAbsorption(group, results);
	CellAvgResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK; 
    }

// Check total opacity results.

    {
	std::cout << " ** Total Opacity **" << std::endl;
	std::cout << " Cell  Opacity (1/m)      Relative Error" 
		  << std::endl;
	std::vector<double> answer(ncell);
	answer[0] = 1.012916666667e+01;
	answer[1] = 5.475416666667e+01;
	answer[2] = 3.038750000000e+01;
	answer[3] = 1.095083333333e+02;
	FTVD results(ncell);
	int group = 1;
	msf.getSigmaTotal(group, results);
	CellAvgResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK; 
    }

// Check emission opacity results.

    {
	std::cout << " ** Emission Opacity **" << std::endl;
	std::cout << " Cell  Opacity (1/m)      Relative Error" 
		  << std::endl;
	std::vector<double> answer(ncell);
	answer[0] = 1.012916666667e+01;
	answer[1] = 5.475416666667e+01;
	answer[2] = 3.038750000000e+01;
	answer[3] = 1.095083333333e+02;
	FTVD results(ncell);
	int group = 1;
	msf.getSigmaEmission(group, results);
	CellAvgResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK; 
    }

// Check electron-ion coupling results.

    {
	std::cout << " ** Electron-Ion Coupling **" << std::endl;
	std::cout << " Cell  eic (J/m**3-s-K)   Relative Error" 
		  << std::endl;
	std::vector<double> answer(ncell);
	answer[0]=   1.342676328483e+14;
	answer[1]=   2.332007110693e+14;
	answer[2]=   4.028028985449e+14;
	answer[3]=   4.664014221386e+14;
	FTVD results(ncell);
	msf.getElectronIonCoupling(results);
	CellAvgResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK; 
    }

// Check electron conduction coefficient results.

    {
	std::cout << " ** Electron Conduction Coefficient **" << std::endl;
	std::cout << " Cell  Cond-Coeff (J/m-s-K) Relative Error" 
		  << std::endl;
	std::vector<double> answer(ncell);
	answer[0]=   1.015572778343e+05;
	answer[1]=   3.007088846946e+05;
	answer[2]=   3.046718335030e+05;
	answer[3]=   6.014177693892e+05;
	FTVD results(ncell);
	msf.getElectronConductionCoeff(results);
	CellAvgResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK; 
    }

// Check ion conduction coefficient results.

    {
	std::cout << " ** Ion Conduction Coefficient **" << std::endl;
	std::cout << " Cell  Cond-Coeff (J/m-s-K) Relative Error" 
		  << std::endl;
	std::vector<double> answer(ncell);
	answer[0]=   5.422036503714e+00;
	answer[1]=   7.900607352896e+00;
	answer[2]=   1.626610951114e+01;
	answer[3]=   1.580121470579e+01;
	FTVD results(ncell);
	msf.getIonConductionCoeff(results);
	CellAvgResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK; 
    }

// Check electron specific heat results.

    {
	std::cout << " ** Electron Specific Heat **" << std::endl;
	std::cout << " Cell  Cve (J/m**3-K)     Relative Error" 
		  << std::endl;
	std::vector<double> answer(ncell);
	answer[0]=   8.358510441380e+00;
	answer[1]=   4.214664257810e+01;
	answer[2]=   2.507553132414e+01;
	answer[3]=   8.429328515620e+01;
	FTVD results(ncell);
	msf.getElectronSpecificHeat(results);
	CellAvgResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK; 
    }

// Check ion specific heat results.

    {
	std::cout << " ** Ion Specific Heat **" << std::endl;
	std::cout << " Cell  Cvi (J/m**3-K)     Relative Error" 
		  << std::endl;
	std::vector<double> answer(ncell);
	answer[0]=   5.950195393002e+00;
	answer[1]=   7.034019568151e+01;
	answer[2]=   1.785058617901e+01;
	answer[3]=   1.406803913630e+02;
	FTVD results(ncell);
	msf.getIonSpecificHeat(results);
	CellAvgResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK;
    }

// Check electron temperature results.

    {
	std::cout << " ** Electron Temperature **" << std::endl;
	std::cout << " Cell  Temperature (K)    Relative Error" 
		  << std::endl;
	std::vector<double> answer(ncell);
	answer[0] = 6.025766197348e+00;
	answer[1] = 1.932154402511e+01;
	answer[2] = 1.807729859204e+01;
    	answer[3] = 3.864308805021e+01;
	FTVD results(ncell);
	msf.getElectronTemperature(results);
	CellAvgResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK; 
    }

// Check ion temperature results.

    {
	std::cout << " ** Ion Temperature **" << std::endl;
	std::cout << " Cell  Temperature (K)    Relative Error" 
		  << std::endl;
	std::vector<double> answer(ncell);
	answer[0]=   5.925879248264e+00;
	answer[1]=   2.171458879996e+01;
	answer[2]=   1.777763774479e+01;
	answer[3]=   4.342917759992e+01;
	FTVD results(ncell);
	msf.getIonTemperature(results);
	CellAvgResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK; 
    }
    
// Check electron specific heat by material results.

    {
	std::cout << " ** Electron Specific Heat By Material **" << std::endl;
	std::cout << " Cell Mat   Cve (J/m**3-K)     Relative Error" 
		  << std::endl;
        FTVVD results(ncell);
        std::vector<std::vector<double> > answer(ncell);
	FTVVD::iterator resit=results.begin();
	for (int icell = 0; icell<ncell; icell++)
	{
	    //	    (*resit++).resize(nmat[icell]);
	    answer[icell].resize(nmat[icell]);
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
	ByMatResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK; 
    }


// Check ion specific heat by material results.

    {
	std::cout << " ** Ion Specific Heat By Material **" << std::endl;
	std::cout << " Cell Mat   Cvi (J/m**3-K)     Relative Error" 
		  << std::endl;
        FTVVD results(ncell);
        std::vector<std::vector<double> > answer(ncell);
	FTVVD::iterator resit=results.begin();
	for (int icell = 0; icell<ncell; icell++)
	{
	    //	    (*resit++).resize(nmat[icell]);
	    answer[icell].resize(nmat[icell]);
	}
	answer[0]  [0]=   7.650146580363e+00;
	answer[0]  [1]=   5.100219799322e+00;
	answer[1]  [0]=   1.530029316073e+01;
	answer[1]  [1]=   1.020043959864e+01;
	answer[1]  [2]=   1.287800005770e+02;
	answer[2]  [0]=   2.295043974109e+01;
	answer[2]  [1]=   1.530065939797e+01;
	answer[3]  [0]=   3.060058632145e+01;
	answer[3]  [1]=   2.040087919729e+01;
	answer[3]  [2]=   2.575600011540e+02;

	msf.getIonSpecificHeatByMat(results);
	ByMatResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK;
    }


// Check electron temperature by material results.

    {
	std::cout << " ** Electron Temperature By Material **" << std::endl;
	std::cout << " Cell Mat   Te (K)             Relative Error" 
		  << std::endl;
        FTVVD results(ncell);
        std::vector<std::vector<double> > answer(ncell);
	FTVVD::iterator resit=results.begin();
	for (int icell = 0; icell<ncell; icell++)
	{
	    // (*resit++).resize(nmat[icell]);
	    answer[icell].resize(nmat[icell]);
	}
	answer[0]  [0]=   3.520000000000e+00;
	answer[0]  [1]=   7.040000000000e+00;
	answer[1]  [0]=   7.040000000000e+00;
	answer[1]  [1]=   1.408000000000e+01;
	answer[1]  [2]=   2.112000000000e+01;
	answer[2]  [0]=   1.056000000000e+01;
	answer[2]  [1]=   2.112000000000e+01;
	answer[3]  [0]=   1.408000000000e+01;
	answer[3]  [1]=   2.816000000000e+01;
	answer[3]  [2]=   4.224000000000e+01;
	msf.getElectronTempByMat(results);
	ByMatResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK;
    }

// Check ion temperature by material results.

    {
	std::cout << " ** Ion Temperature By Material **" << std::endl;
	std::cout << " Cell Mat   Ti (K)             Relative Error" 
		  << std::endl;
        FTVVD results(ncell);
        std::vector<std::vector<double> > answer(ncell);
	FTVVD::iterator resit=results.begin();
	for (int icell = 0; icell<ncell; icell++)
	{
	    //	    (*resit++).resize(nmat[icell]);
	    answer[icell].resize(nmat[icell]);
	}
	answer[0]  [0]=   3.771000000000e+00;
	answer[0]  [1]=   7.542000000000e+00;
	answer[1]  [0]=   7.542000000000e+00;
	answer[1]  [1]=   1.508400000000e+01;
	answer[1]  [2]=   2.262600000000e+01;
	answer[2]  [0]=   1.131300000000e+01;
	answer[2]  [1]=   2.262600000000e+01;
	answer[3]  [0]=   1.508400000000e+01;
	answer[3]  [1]=   3.016800000000e+01;
	answer[3]  [2]=   4.525200000000e+01;
	msf.getIonTempByMat(results);
	ByMatResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK;
    }

    // Check Density results

    {
	std::cout << " ** Density By Material **" << std::endl;
	std::cout << " Cell Mat   Dens (kg/m**3))    Relative Error" 
		  << std::endl;
        FTVVD results(ncell);
        std::vector<std::vector<double> > answer(ncell);
	FTVVD::iterator resit=results.begin();
	for (int icell = 0; icell<ncell; icell++)
	{
	    //	    (*resit++).resize(nmat[icell]);
	    answer[icell].resize(nmat[icell]);
	}
	answer[0]  [0]=   2.125000000000e+00;
	answer[0]  [1]=   4.250000000000e+00;
	answer[1]  [0]=   4.250000000000e+00;
	answer[1]  [1]=   8.500000000000e+00;
	answer[1]  [2]=   1.275000000000e+01;
	answer[2]  [0]=   6.375000000000e+00;
	answer[2]  [1]=   1.275000000000e+01;
	answer[3]  [0]=   8.500000000000e+00;
	answer[3]  [1]=   1.700000000000e+01;
	answer[3]  [2]=   2.550000000000e+01;
	msf.getDensity(results);
	ByMatResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK;
    }

    // Check Volume Fraction results

    {
	std::cout << " ** Material Volume Fractions **" << std::endl;
	std::cout << " Cell Mat   Volume Frac (nd)   Relative Error" 
		  << std::endl;
        FTVVD results(ncell);
        std::vector<std::vector<double> > answer(ncell);
	FTVVD::iterator resit=results.begin();
	for (int icell = 0; icell<ncell; icell++)
	{
	    //	    (*resit++).resize(nmat[icell]);
	    answer[icell].resize(nmat[icell]);
	}
	answer[0]  [0]=   3.333333333333e-01;
	answer[0]  [1]=   6.666666666667e-01;
	answer[1]  [0]=   1.666666666667e-01;
	answer[1]  [1]=   3.333333333333e-01;
	answer[1]  [2]=   5.000000000000e-01;
	answer[2]  [0]=   3.333333333333e-01;
	answer[2]  [1]=   6.666666666667e-01;
	answer[3]  [0]=   1.666666666667e-01;
	answer[3]  [1]=   3.333333333333e-01;
	answer[3]  [2]=   5.000000000000e-01;
	msf.getVolumeFraction(results);
	ByMatResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK;
    }

    // Check Material ID results

    {
	std::cout << " ** Material ID's **" << std::endl;
	std::cout << " Cell Mat  ID Error" 
		  << std::endl;
        FTVVI results(ncell);
        std::vector<std::vector<int> > answer(ncell);
	FTVVI::iterator resit=results.begin();
	for (int icell = 0; icell<ncell; icell++)
	{
	    //	    (*resit++).resize(nmat[icell]);
	    answer[icell].resize(nmat[icell]);
	}
	answer[0]  [0]=   1;
	answer[0]  [1]=   2;
	answer[1]  [0]=   1;
	answer[1]  [1]=   2;
	answer[1]  [2]=   3;
	answer[2]  [0]=   1;
	answer[2]  [1]=   2;
	answer[3]  [0]=   1;
	answer[3]  [1]=   2;
	answer[3]  [2]=   3;
	msf.getMatId(results);
	ByMatResultsOK(results, answer, eps, results_OK);
	passed = passed && results_OK;
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
