//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/test/tstAnalytic_Gray_Opacity.cc
 * \author Thomas M. Evans
 * \date   Mon Sep 24 12:08:55 2001
 * \brief  Analytic_Gray_Opacity test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_analytic_test.hh"
#include "../Release.hh"
#include "../Analytic_Gray_Opacity.hh"
#include "../Analytic_Models.hh"
#include "cdi/CDI.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <vector>
#include <iostream>
#include <string>
#include <typeinfo>
#include <algorithm>
#include <cmath>
#include <sstream>

using namespace std;

using rtt_cdi_analytic::Analytic_Gray_Opacity;
using rtt_cdi_analytic::Analytic_Opacity_Model;
using rtt_cdi_analytic::Constant_Analytic_Opacity_Model;
using rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model;
using rtt_cdi::CDI;
using rtt_cdi::GrayOpacity;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
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
    if (grayp->getModelType() != rtt_cdi::ANALYTIC)          ITFAILS;
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

    vector<double> opacities = grayp->getOpacity(T, rho[0]);
    if (opacities.size() != 6) ITFAILS;

    for (int i = 0; i < T.size(); i++)
    {
	double ref         = 10.0 / (T[i]*T[i]*T[i]);
	double error       = fabs(grayp->getOpacity(T[i], rho[i]) - ref);
	double error_field = fabs(opacities[i] - ref);  

	if (error > 1.0e-12 * ref)       ITFAILS; 
	if (error_field > 1.0e-12 * ref) ITFAILS; 
    }

    // check to make sure we can't unpack an unregistered analytic model
    vector<char> packed = anal_opacity.pack();

    bool caught = false;
    try
    {
	Analytic_Gray_Opacity ngray(packed);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	caught = true;
	ostringstream message;
	message << "Caught the following assertion: " << ass.what();
	PASSMSG(message.str());
    }
    if (!caught)
    {
	FAILMSG("Failed to catch unregistered analyic model assertion");
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

    absorption = new const Analytic_Gray_Opacity(amodel, rtt_cdi::ABSORPTION);
    scattering = new const Analytic_Gray_Opacity(smodel, rtt_cdi::SCATTERING);

    if (!absorption) FAILMSG("Failed to build absorption analytic opacity")
    if (!scattering) FAILMSG("Failed to build scattering analytic opacity")


    // make a CDI for scattering and absorption
    CDI cdi;
    cdi.setGrayOpacity(scattering);
    cdi.setGrayOpacity(absorption);

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
	double ref = 100.0 / (T[i]*T[i]*T[i]);
	rtt_cdi::Model model   = rtt_cdi::ANALYTIC;
	rtt_cdi::Reaction abs  = rtt_cdi::ABSORPTION;
	rtt_cdi::Reaction scat = rtt_cdi::SCATTERING; 

	double error = fabs(cdi.gray(model,abs)->getOpacity(T[i], rho[i])
			    - ref);

	if (error > 1.0e-12 * ref) ITFAILS; 

	error = fabs(cdi.gray(model, scat)->getOpacity(T[i], rho[i]) - 1.0);

	if (error > 1.0e-12)       ITFAILS;
    }
}

//---------------------------------------------------------------------------//

void packing_test()
{
    // test the packing
    vector<char> packed;
    {
	// lets make two models
	SP<Analytic_Opacity_Model> amodel
	    (new Polynomial_Analytic_Opacity_Model(0.0,100.0,-3.0,1.0,0.0));
	
	Analytic_Gray_Opacity absorption(amodel, rtt_cdi::ABSORPTION);

	packed = absorption.pack();
    }

    // now unpack and test
    Analytic_Gray_Opacity ngray(packed);

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
	double ref = 100.0 / (T[i]*T[i]*T[i]);

	double error = fabs(ngray.getOpacity(T[i], rho[i]) - ref);

	if (error > 1.0e-12 * ref) ITFAILS; 
    }
    
    if (ngray.getReactionType() != rtt_cdi::ABSORPTION) ITFAILS;
    if (ngray.getModelType() != rtt_cdi::ANALYTIC)      ITFAILS;

    if (rtt_cdi_analytic_test::passed)
	PASSMSG("Analytic_Gray_Opacity packing test passes.");
}

//---------------------------------------------------------------------------//

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

    try
    {
	// >>> UNIT TESTS
	constant_test();
	user_defined_test();
	CDI_test();

	packing_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstAnalytic_Gray_Opacity, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_cdi_analytic_test::passed) 
    {
        cout << "**** tstAnalytic_Gray_Opacity Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstAnalytic_Gray_Opacity." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstAnalytic_Gray_Opacity.cc
//---------------------------------------------------------------------------//
