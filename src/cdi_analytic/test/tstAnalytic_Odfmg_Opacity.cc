//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/test/tstAnalytic_Odfmg_Opacity.cc
 * \author Thomas M. Evans
 * \date   Tue Nov 13 17:24:12 2001
 * \brief  Analytic_Odfmg_Opacity test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_analytic_test.hh"
#include "ds++/Release.hh"
#include "../Analytic_Odfmg_Opacity.hh"
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

using rtt_cdi_analytic::Analytic_Odfmg_Opacity;
using rtt_cdi_analytic::Analytic_Opacity_Model;
using rtt_cdi_analytic::Constant_Analytic_Opacity_Model;
using rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model;
using rtt_cdi::CDI;
using rtt_cdi::OdfmgOpacity;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

bool checkOpacityEquivalence( vector< vector<double> > sigma,
                              vector<double> ref)
{
    bool itPasses = true;

    if (sigma.size() != ref.size())
    {
        cout << "Mismatch in number of groups: reference "
             << ref.size() << ", from opacity " << sigma.size() << endl;

        return false;
    }

    for (size_t group = 0; group < sigma.size(); group++)
    {
        for (size_t band = 0; band < sigma[group].size(); band++)
        {
            if (!soft_equiv(sigma[group][band],ref[group]))
            {
                itPasses = false;
                cout << "Mismatch in opacities for group " << group
                     << " band " << band << endl;
            }
        }
    }
    return itPasses;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void odfmg_test()
{
    // group structure
    vector<double> groups(4, 0.0);
    groups[0] = 0.05;
    groups[1] = 0.5;
    groups[2] = 5.0;
    groups[3] = 50.0;

    // band strucutre
    vector<double> bands(3, 0.0);
    bands[0] = 0.0;
    bands[1] = 0.75;
    bands[2] = 1.0;

    vector<SP<Analytic_Opacity_Model> > models(3);

    // make a Marshak (user-defined) model for the first group
    models[0] = new rtt_cdi_analytic_test::Marshak_Model(100.0);

    // make a Polynomial model for the second group
    models[1] = new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
        1.5, 0.0, 0.0, 0.0);

    // make a Constant model for the third group
    models[2] = new rtt_cdi_analytic::Constant_Analytic_Opacity_Model(3.0);

    // make an analytic multigroup opacity object for absorption
    Analytic_Odfmg_Opacity opacity(groups, bands, models, rtt_cdi::ABSORPTION);

    // check the interface to multigroup opacity
    {
        string desc = "Analytic Odfmg Absorption";

        if (opacity.getOpacityModelType()!=
            rtt_cdi::ANALYTIC_TYPE)                           ITFAILS;
        if (opacity.data_in_tabular_form())                   ITFAILS;
        if (opacity.getReactionType() != rtt_cdi::ABSORPTION) ITFAILS;
        if (opacity.getModelType() != rtt_cdi::ANALYTIC)      ITFAILS;
        if (opacity.getNumTemperatures() != 0)                ITFAILS;
        if (opacity.getNumDensities() != 0)                   ITFAILS;
        if (opacity.getTemperatureGrid() != vector<double>()) ITFAILS;
        if (opacity.getDensityGrid() != vector<double>())     ITFAILS;
        if (opacity.getNumGroups() != 3)                      ITFAILS;
        if (opacity.getNumGroupBoundaries() != 4)             ITFAILS;
        if (opacity.getNumBands() != 2)                       ITFAILS;
        if (opacity.getNumBandBoundaries() != 3)              ITFAILS;
        if (opacity.getEnergyPolicyDescriptor() != "odfmg")   ITFAILS;
        if (opacity.getDataDescriptor() != desc)              ITFAILS;
        if (opacity.getDataFilename() != string())            ITFAILS;

        if (opacity.get_Analytic_Model(1, 1) != models[0])    ITFAILS;
        if (opacity.get_Analytic_Model(2, 1) != models[1])    ITFAILS;
        if (opacity.get_Analytic_Model(3, 1) != models[2])    ITFAILS;
    }

    {
        // make an analytic multigroup opacity object for scattering
        Analytic_Odfmg_Opacity opacity(groups,
                                       bands,
                                       models,
                                       rtt_cdi::SCATTERING);
        string desc = "Analytic Odfmg Scattering";

        if (opacity.getDataDescriptor() != desc)              ITFAILS;
    }
    {
        // make an analytic multigroup opacity object for scattering
        Analytic_Odfmg_Opacity opacity(groups, bands, models, rtt_cdi::TOTAL);
        string desc = "Analytic Odfmg Total";

        if (opacity.getDataDescriptor() != desc)              ITFAILS;
    }

    // check the group structure
    vector<double> get_groups = opacity.getGroupBoundaries();

    if (soft_equiv(get_groups.begin(), get_groups.end(), 
                   groups.begin(), groups.end()))
    {
        PASSMSG("Group boundaries match.");
    }
    else
    {
        FAILMSG("Group boundaries do not match.");
    }

    // check the band structure
    vector<double> get_bands = opacity.getBandBoundaries();

    if (soft_equiv(get_bands.begin(), get_bands.end(), 
                   bands.begin(), bands.end()))
    {
        PASSMSG("Band boundaries match.");
    }
    else
    {
        FAILMSG("Band boundaries do not match.");
    }


    // >>> get opacities

    // scalar density and temperature
    vector<double> ref(3, 0.0);
    ref[0] = 100.0 / 8.0;
    ref[1] = 1.5;
    ref[2] = 3.0;
        
    // load groups * bands opacities; all bands inside each group should be
    // the same
    vector< vector<double> > sigma = opacity.getOpacity(2.0, 3.0); 

    // check for each band and group
    bool itPasses = checkOpacityEquivalence(sigma, ref);

    if (itPasses)
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
    vector< vector< vector<double> > > sig_t
        = opacity.getOpacity(data_field, 3.0);
    vector< vector< vector<double> > > sig_rho
        = opacity.getOpacity(2.0, data_field);

    itPasses = true;
    for (int i = 0; i < 3; i++)
    {
        itPasses = itPasses && checkOpacityEquivalence(sig_t[i], ref);
    }

    if (itPasses)
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

    itPasses = true;
    for (int i = 0; i < 3; i++)
    {
        itPasses = itPasses && checkOpacityEquivalence(sig_rho[i], ref);
    }

    if (itPasses)
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

    // Test the get_Analytic_Model() member function.
    SP<Analytic_Opacity_Model const> my_mg_opacity_model
        = opacity.get_Analytic_Model(1);
    SP<Analytic_Opacity_Model const> expected_model( models[0] );

    if( expected_model == my_mg_opacity_model )
        PASSMSG("get_Analytic_Model() returned the expected MG Opacity model.")
    else
        FAILMSG("get_Analytic_Model() did not return the expected MG Opacity model.")
    
    return;
}

//---------------------------------------------------------------------------//

void test_CDI()
{
    // group structure
    vector<double> groups(4, 0.0);
    groups[0] = 0.05;
    groups[1] = 0.5;
    groups[2] = 5.0;
    groups[3] = 50.0;

    // band strucutre
    vector<double> bands(3, 0.0);
    bands[0] = 0.0;
    bands[1] = 0.75;
    bands[2] = 1.0;

    vector<SP<Analytic_Opacity_Model> > models(3);

    // make a Marshak (user-defined) model for the first group
    models[0] = new rtt_cdi_analytic_test::Marshak_Model(100.0);

    // make a Polynomial model for the second group
    models[1] = new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
        1.5, 0.0, 0.0, 0.0);

    // make a Constant model for the third group
    models[2] = new rtt_cdi_analytic::Constant_Analytic_Opacity_Model(3.0);

    // make an analytic multigroup opacity object for absorption
    SP<const OdfmgOpacity> odfmg( 
        new Analytic_Odfmg_Opacity(groups, bands, models, rtt_cdi::ABSORPTION));

    // make a CDI object
    CDI cdi;

    // set the multigroup opacity
    cdi.setOdfmgOpacity(odfmg);

    // check the energy groups from CDI
    vector<double> odfmg_groups = CDI::getFrequencyGroupBoundaries();

    if (soft_equiv(odfmg_groups.begin(), odfmg_groups.end(), 
                   groups.begin(), groups.end()))
    {
        PASSMSG("CDI Group boundaries match.");
    }
    else
    {
        FAILMSG("CDI Group boundaries do not match.");
    }
        
    // check the energy groups from CDI
    vector<double> odfmg_bands = CDI::getOpacityCdfBandBoundaries();

    if (soft_equiv(odfmg_bands.begin(), odfmg_bands.end(), 
                   bands.begin(), bands.end()))
    {
        PASSMSG("CDI band boundaries match.");
    }
    else
    {
        FAILMSG("CDI band boundaries do not match.");
    }


    // do a quick access test for getOpacity

    // scalar density and temperature
    vector< vector<double> > sigma
        = cdi.odfmg(rtt_cdi::ANALYTIC, 
                    rtt_cdi::ABSORPTION)->getOpacity(2.0, 3.0);

    vector<double> ref(3, 0.0);
    ref[0] = 100.0 / 8.0;
    ref[1] = 1.5;
    ref[2] = 3.0;

    if (checkOpacityEquivalence(sigma, ref))
    {
        ostringstream message;
        message << "CDI odfmg opacities are correct for "
                << "scalar temperature and scalar density.";
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message << "CDI odfmg opacities are NOT correct for "
                << "scalar temperature and scalar density.";
        FAILMSG(message.str());
    }
}

//---------------------------------------------------------------------------//

void packing_test()
{
    vector<char> packed;

    // group structure
    vector<double> groups(4, 0.0);
    groups[0] = 0.05;
    groups[1] = 0.5;
    groups[2] = 5.0;
    groups[3] = 50.0;
        
    // band strucutre
    vector<double> bands(3, 0.0);
    bands[0] = 0.0;
    bands[1] = 0.75;
    bands[2] = 1.0;

    {
        vector<SP<Analytic_Opacity_Model> > models(3);

        // make a Polynomial model for the first group
        models[0] = new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
            0.0, 100.0, -3.0, 0.0);

        // make a Polynomial model for the second group
        models[1] = new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
            1.5, 0.0, 0.0, 0.0);

        // make a Constant model for the third group
        models[2] = new rtt_cdi_analytic::Constant_Analytic_Opacity_Model(3.0);

        // make an analytic multigroup opacity object for absorption
        SP<const OdfmgOpacity> odfmg
            (new Analytic_Odfmg_Opacity(groups, bands, models, 
                                        rtt_cdi::ABSORPTION));

        // pack it
        packed = odfmg->pack();
    }

    // now unpack it
    Analytic_Odfmg_Opacity opacity(packed);

    // now check it

    // check the interface to multigroup opacity
    {
        string desc = "Analytic Odfmg Absorption";

        if (opacity.data_in_tabular_form())                   ITFAILS;
        if (opacity.getReactionType() != rtt_cdi::ABSORPTION) ITFAILS;
        if (opacity.getModelType() != rtt_cdi::ANALYTIC)      ITFAILS;
        if (opacity.getNumTemperatures() != 0)                ITFAILS;
        if (opacity.getNumDensities() != 0)                   ITFAILS;
        if (opacity.getTemperatureGrid() != vector<double>()) ITFAILS;
        if (opacity.getDensityGrid() != vector<double>())     ITFAILS;
        if (opacity.getNumGroups() != 3)                      ITFAILS;
        if (opacity.getNumGroupBoundaries() != 4)             ITFAILS;
        if (opacity.getEnergyPolicyDescriptor() != "odfmg")   ITFAILS;
        if (opacity.getDataDescriptor() != desc)              ITFAILS;
        if (opacity.getDataFilename() != string())            ITFAILS;
    }

    // check the group structure
    vector<double> mg_groups = opacity.getGroupBoundaries();

    if (soft_equiv(mg_groups.begin(), mg_groups.end(), 
                   groups.begin(), groups.end()))
    {
        PASSMSG("Group boundaries for unpacked ODFMG opacity match.");
    }
    else
    {
        FAILMSG("Group boundaries for unpacked ODFMG do not match.");
    }
        
    // check the band structure
    vector<double> get_bands = opacity.getBandBoundaries();

    if (soft_equiv(get_bands.begin(), get_bands.end(), 
                   bands.begin(), bands.end()))
    {
        PASSMSG("Band boundaries match.");
    }
    else
    {
        FAILMSG("Band boundaries do not match.");
    }

    // >>> get opacities

    // scalar density and temperature
    vector< vector<double> > sigma = opacity.getOpacity(2.0, 3.0);
    vector<double> ref(3, 0.0);
    ref[0] = 100.0 / 8.0;
    ref[1] = 1.5;
    ref[2] = 3.0;

    if (checkOpacityEquivalence(sigma, ref))
    {
        ostringstream message;
        message << "Analytic multigroup opacities for unpacked MG opacity "
                << "are correct for "
                << "scalar temperature and scalar density.";
        PASSMSG(message.str());
    }
    else
    {
        ostringstream message;
        message << "Analytic multigroup opacities for unpacked MG opacity "
                << "are NOT correct for "
                << "scalar temperature and scalar density.";
        FAILMSG(message.str());
    }

    // make sure we catch an assertion showing that we cannot unpack an
    // unregistered opacity
    {
        vector<SP<Analytic_Opacity_Model> > models(3);

        // make a Marshak (user-defined) model for the first group
        models[0] = new rtt_cdi_analytic_test::Marshak_Model(100.0);

        // make a Polynomial model for the second group
        models[1] = new rtt_cdi_analytic::Polynomial_Analytic_Opacity_Model(
            1.5, 0.0, 0.0, 0.0);

        // make a Constant model for the third group
        models[2] = new rtt_cdi_analytic::Constant_Analytic_Opacity_Model(3.0);

        // make an analytic multigroup opacity object for absorption
        SP<const OdfmgOpacity> odfmg( 
            new Analytic_Odfmg_Opacity(groups, bands, models,
                                       rtt_cdi::ABSORPTION));

        packed = odfmg->pack();
    }

    // we should catch an assertion when unpacking this because the
    // Marshak_Model is not registered in rtt_cdi::Opacity_Models
    bool caught = false;
    try
    {
        Analytic_Odfmg_Opacity nmg(packed);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
        caught = true;
        ostringstream message;
        message << "Caught the following assertion, " << ass.what();
        PASSMSG(message.str());
    }
    if (!caught)
    {
        FAILMSG("Failed to catch unregistered analyic model assertion");
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (string(argv[arg]) == "--version")
        {
            cout << argv[0] << ": version " << rtt_dsxx::release() 
                 << endl;
            return 0;
        }

    try
    {
        // >>> UNIT TESTS
        odfmg_test();

        test_CDI();

        packing_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
        cout << "While testing tstAnalytic_Odfmg_Opacity, " << ass.what()
             << endl;
        return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_cdi_analytic_test::passed) 
    {
        cout << "**** tstAnalytic_Odfmg_Opacity Test: PASSED" 
             << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstAnalytic_Odfmg_Opacity." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstAnalytic_Odfmg_Opacity.cc
//---------------------------------------------------------------------------//
