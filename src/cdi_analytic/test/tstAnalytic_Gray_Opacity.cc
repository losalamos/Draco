//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/test/tstAnalytic_Gray_Opacity.cc
 * \author Thomas M. Evans
 * \date   Wed Aug 29 14:38:57 2001
 * \brief  Analytic_Gray_Opacity test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "CDI_Test.hh"
#include "../Release.hh"
#include "../Analytic_Gray_Opacity.hh"
#include "../Analytic_Models.hh"
#include "cdi/CDI.hh"
#include "cdi/OpacityCommon.hh"
#include "cdi/GrayOpacity.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <vector>
#include <iostream>
#include <string>
#include <typeinfo>
#include <algorithm>
#include <cmath>

using namespace std;

using rtt_cdi_analytic::Analytic_Gray_Opacity;
using rtt_cdi_analytic::Analytic_Opacity_Model;
using rtt_cdi_analytic::Constant_Analytic_Opacity_Model;
using rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model;
using rtt_cdi::CDI;
using rtt_cdi::GrayOpacity;
using rtt_dsxx::SP;

bool passed = true;
#define ITFAILS passed = rtt_cdi_analytic_test::fail(__LINE__);

//---------------------------------------------------------------------------//

void constant_test()
{
    // make an analytic gray opacity that returns the total opacity for a
    // constant model
    const double constant_opacity = 5.0;

    SP<Analytic_Opacity_Model> model
	(new Constant_Analytic_Opacity_Model(constant_opacity));

    Analytic_Gray_Opacity anal_opacity(model, rtt_cdi::TOTAL);

    // link to GrayOpacity
    GrayOpacity *grayp = &anal_opacity;

    // Some basic tests
    if (grayp->data_in_tabular_form())                       ITFAILS;
    if (grayp->getEnergyPolicyDescriptor() != "gray")        ITFAILS;
    if (grayp->getDataDescriptor() != "Analytic Gray Total") ITFAILS;
    if (grayp->getReactionType() != rtt_cdi::TOTAL)          ITFAILS;
    if (typeid(grayp) != typeid(GrayOpacity *))              ITFAILS;
    if (typeid(*grayp) != typeid(Analytic_Gray_Opacity))     ITFAILS;

    // check the output
    vector<double> T(10);
    vector<double> rho(10);

    for (int i = 0; i < T.size(); i++)
    {
	T[i]   = 0.1 + i/100.0;
	rho[i] = 1.0 + i/10.0;

	if (grayp->getOpacity(T[i], rho[i]) != constant_opacity) ITFAILS;
    }

    vector<double> opacity_T   = grayp->getOpacity(T, 3.0);
    vector<double> opacity_rho = grayp->getOpacity(1.0, rho);
    vector<double> ref(10, constant_opacity);
    
    if (opacity_T != ref)   ITFAILS;
    if (opacity_rho != ref) ITFAILS;
}

//---------------------------------------------------------------------------//

void user_defined_test()
{
    // make the user defined Marshak model
    SP<Analytic_Opacity_Model> model
	(new rtt_cdi_analytic_test::Marshak_Model(10.0));

    Analytic_Gray_Opacity anal_opacity(model, rtt_cdi::TOTAL);
    GrayOpacity *grayp = &anal_opacity;

    vector<double> T(6);
    vector<double> rho(6);
    {
	T[0] = .993;
	T[1] = .882;
	T[2] = .590;
	T[3] = .112;
	T[4] = .051;
	T[5] = .001;

	std::fill(rho.begin(), rho.end(), 3.0);
    }

    for (int i = 0; i < T.size(); i++)
    {
	double ref   = 10.0 / (T[i]*T[i]*T[i]);
	double error = fabs(grayp->getOpacity(T[i], rho[i]) - ref);

	if (error > 1.0e-12 * ref) ITFAILS; 
    }
}

//---------------------------------------------------------------------------//

void CDI_test()
{
    // lets make a marshak model gray opacity for scattering and absorption
    SP<const GrayOpacity> absorption;
    SP<const GrayOpacity> scattering;

    // lets make two models
    SP<Analytic_Opacity_Model> amodel
	(new Polynomial_Analytic_Opacity_Model(0.0,100.0,-3.0,1.0,0.0));
    SP<Analytic_Opacity_Model> smodel
	(new Constant_Analytic_Opacity_Model(1.0));

    absorption = new Analytic_Gray_Opacity(amodel, rtt_cdi::ABSORPTION);
    scattering = new Analytic_Gray_Opacity(smodel, rtt_cdi::SCATTERING);

    // make a CDI for scattering and absorption (should you be able to add
    // multiple opacities to a CDI)
    CDI scat_CDI(scattering);
    CDI abs_CDI(absorption);

    // now check some data
    vector<double> T(6);
    vector<double> rho(6, 3.0);
    {
	T[0] = .993;
	T[1] = .882;
	T[2] = .590;
	T[3] = .112;
	T[4] = .051;
	T[5] = .001;

	std::fill(rho.begin(), rho.end(), 3.0);
    }

    for (int i = 0; i < T.size(); i++)
    {
	double ref   = 100.0 / (T[i]*T[i]*T[i]);
	double error = fabs(abs_CDI.gray()->getOpacity(T[i], rho[i]) - ref);

	if (error > 1.0e-12 * ref) ITFAILS; 
    }
    
}

//---------------------------------------------------------------------------//
// main

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_cdi_analytic::release()
		 << endl;
	    return 0;
	}

    // some simple tests
    try
    {
	constant_test();
	user_defined_test();
	CDI_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "Caught: " << ass.what() << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "****************************************************" 
	 << endl;
    if (passed) 
    {
        cout << "**** tstAnalytic_Gray_Opacity Self Test: PASSED ****" 
	     << endl;
    }
    cout <<     "****************************************************" 
	 << endl;
    cout << endl;

    cout << "Done testing tstAnalytic_Gray_Opacity." << endl;
}

//---------------------------------------------------------------------------//
//                              end of tstAnalytic_Gray_Opacity.cc
//---------------------------------------------------------------------------//
