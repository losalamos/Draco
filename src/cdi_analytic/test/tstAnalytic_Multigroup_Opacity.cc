//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/test/tstAnalytic_Multigroup_Opacity.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 13 17:24:12 2001
 * \brief  Analytic_Multigroup_Opacity test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_analytic_test.hh"
#include "../Release.hh"
#include "../Analytic_Multigroup_Opacity.hh"
#include "../Analytic_Models.hh"
#include "cdi/CDI.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <typeinfo>
#include <algorithm>
#include <sstream>

using namespace std;

using rtt_cdi_analytic::Analytic_Multigroup_Opacity;
using rtt_cdi_analytic::Analytic_Opacity_Model;
using rtt_cdi_analytic::Constant_Analytic_Opacity_Model;
using rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model;
using rtt_cdi::CDI;
using rtt_cdi::MultigroupOpacity;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void multigroup_test()
{
    // group structure
    vector<double> groups(4, 0.0);
    {
	groups[0] = 0.05;
	groups[1] = 0.5;
	groups[2] = 5.0;
	groups[3] = 50.0;
    }
    
    vector<SP<Analytic_Opacity_Model> > models(3);

    // make a Marshak (user-defined) model for the first group
    models[0] = new rtt_cdi_analytic_test::Marshak_Model(100.0);

    // make a Polynomial model for the second group
    models[1] = new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
	1.5, 0.0, 0.0, 1.0, 0.0);

    // make a Constant model for the third group
    models[2] = new rtt_cdi_analytic::Constant_Analytic_Opacity_Model(3.0);

    // make an analytic multigroup opacity object for absorption
    Analytic_Multigroup_Opacity opacity(groups, models, rtt_cdi::ABSORPTION);

    // check the interface to multigroup opacity
    {
	string desc = "Analytic Multigroup Absorption";

	if (opacity.data_in_tabular_form())                   ITFAILS;
	if (opacity.getReactionType() != rtt_cdi::ABSORPTION) ITFAILS;
	if (opacity.getModelType() != rtt_cdi::ANALYTIC)      ITFAILS;
	if (opacity.getNumTemperatures() != 0)                ITFAILS;
	if (opacity.getNumDensities() != 0)                   ITFAILS;
	if (opacity.getTemperatureGrid() != vector<double>()) ITFAILS;
	if (opacity.getDensityGrid() != vector<double>())     ITFAILS;
	if (opacity.getNumGroups() != 3)                      ITFAILS;
	if (opacity.getNumGroupBoundaries() != 4)             ITFAILS;
	if (opacity.getEnergyPolicyDescriptor() != "mg")      ITFAILS;
	if (opacity.getDataDescriptor() != desc)              ITFAILS;
	if (opacity.getDataFilename() != string())            ITFAILS;
    }

    // check the group structure
    vector<double> mg_groups = opacity.getGroupBoundaries();

    if (soft_equiv(mg_groups.begin(), mg_groups.end(), 
		   groups.begin(), groups.end()))
    {
	PASSMSG("Group boundaries match.");
    }
    else
    {
	FAILMSG("Group boundaries do not match.");
    }

    // >>> get opacities
    
    // scalar density and temperature
    vector<double> sigma = opacity.getOpacity(2.0, 3.0);
    vector<double> ref(3, 0.0);
    {
	ref[0] = 100.0 / 8.0;
	ref[1] = 1.5;
	ref[2] = 3.0;
    }
    if (soft_equiv(sigma.begin(), sigma.end(), ref.begin(), ref.end()))
    {
	ostringstream message;
	message << "Analytic multigroup opacities are correct for "
		<< "scalar temperature and scalar density.";
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Analytic multigroup opacities are NOT correct for "
		<< "scalar temperature and scalar density.";
	FAILMSG(message.str());
    }

    // scalar density/temperature + vector density/temperature
    vector<double>          data_field(3, 2.0);
    vector<vector<double> > sig_t   = opacity.getOpacity(data_field, 3.0);
    vector<vector<double> > sig_rho = opacity.getOpacity(2.0, data_field);

    for (int i = 0; i < 3; i++)
    {
	vector<double> &test = sig_t[i];

	if (soft_equiv(test.begin(), test.end(), 
		       ref.begin(), ref.end()))
	{
	    ostringstream message;
	    message << "Analytic multigroup opacities are correct for "
		    << "temperature field and scalar density.";
	    PASSMSG(message.str());
	}
	else
	{
	    ostringstream message;
	    message << "Analytic multigroup opacities are NOT correct for "
		    << "temperature field and scalar density.";
	    FAILMSG(message.str());
	}

	test = sig_rho[i];

	if (soft_equiv(test.begin(), test.end(), 
		       ref.begin(), ref.end()))
	{
	    ostringstream message;
	    message << "Analytic multigroup opacities are correct for "
		    << "density field and scalar temperature.";
	    PASSMSG(message.str());
	}
	else
	{
	    ostringstream message;
	    message << "Analytic multigroup opacities are NOT correct for "
		    << "density field and scalar temperature.";
	    FAILMSG(message.str());
	}
    }
}

//---------------------------------------------------------------------------//

void test_CDI()
{
    // group structure
    vector<double> groups(4, 0.0);
    {
	groups[0] = 0.05;
	groups[1] = 0.5;
	groups[2] = 5.0;
	groups[3] = 50.0;
    }
    
    vector<SP<Analytic_Opacity_Model> > models(3);

    // make a Marshak (user-defined) model for the first group
    models[0] = new rtt_cdi_analytic_test::Marshak_Model(100.0);

    // make a Polynomial model for the second group
    models[1] = new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
	1.5, 0.0, 0.0, 1.0, 0.0);

    // make a Constant model for the third group
    models[2] = new rtt_cdi_analytic::Constant_Analytic_Opacity_Model(3.0);

    // make an analytic multigroup opacity object for absorption
    SP<const MultigroupOpacity> mg( 
	new Analytic_Multigroup_Opacity(groups, models, rtt_cdi::ABSORPTION));

    // make a CDI object
    CDI cdi;

    // set the multigroup opacity
    cdi.setMultigroupOpacity(mg);

    // check the energy groups from CDI
    vector<double> mg_groups = CDI::getFrequencyGroupBoundaries();

    if (soft_equiv(mg_groups.begin(), mg_groups.end(), 
		   groups.begin(), groups.end()))
    {
	PASSMSG("CDI Group boundaries match.");
    }
    else
    {
	FAILMSG("CDI Group boundaries do not match.");
    }

    // do a quick access test for getOpacity
    
    // scalar density and temperature
    vector<double> sigma = cdi.mg(rtt_cdi::ANALYTIC, 
				  rtt_cdi::ABSORPTION)->getOpacity(2.0, 3.0);
    vector<double> ref(3, 0.0);
    {
	ref[0] = 100.0 / 8.0;
	ref[1] = 1.5;
	ref[2] = 3.0;
    }
    if (soft_equiv(sigma.begin(), sigma.end(), ref.begin(), ref.end()))
    {
	ostringstream message;
	message << "CDI multigroup opacities are correct for "
		<< "scalar temperature and scalar density.";
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "CDI multigroup opacities are NOT correct for "
		<< "scalar temperature and scalar density.";
	FAILMSG(message.str());
    }
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
	multigroup_test();

	test_CDI();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstAnalytic_Multigroup_Opacity, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_cdi_analytic_test::passed) 
    {
        cout << "**** tstAnalytic_Multigroup_Opacity Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstAnalytic_Multigroup_Opacity." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstAnalytic_Multigroup_Opacity.cc
//---------------------------------------------------------------------------//
