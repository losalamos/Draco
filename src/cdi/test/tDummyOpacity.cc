//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/tDummyOpacity.cc
 * \author Thomas M. Evans
 * \date   Tue Oct  9 15:50:53 2001
 * \brief  GrayOpacity and Multigroup opacity test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_test.hh"
#include "DummyGrayOpacity.hh"
#include "DummyMultigroupOpacity.hh"
#include "../Release.hh"
#include "../GrayOpacity.hh"
#include "../MultigroupOpacity.hh"
#include "../OpacityCommon.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

using rtt_cdi_test::match;
using rtt_cdi::GrayOpacity;
using rtt_cdi::MultigroupOpacity;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void simple_tests()
{
    // make SPs to gray and multigroup opacities
    SP<GrayOpacity>       gray;
    SP<MultigroupOpacity> mg;

    // Assign and check gray opacity  
    SP<rtt_cdi_test::DummyGrayOpacity> gray_total;
    SP<rtt_cdi_test::DummyGrayOpacity> gray_abs;
    gray_total = new rtt_cdi_test::DummyGrayOpacity();
    gray_abs   = new rtt_cdi_test::DummyGrayOpacity(rtt_cdi::ABSORPTION);

    // check gray opacity for total opacities
    {
	gray = gray_total;
	
	if (!gray->data_in_tabular_form())             ITFAILS;
	if (gray->getReactionType() != rtt_cdi::TOTAL) ITFAILS;
	if (gray->getModelType() != rtt_cdi::ANALYTIC) ITFAILS;
	
	vector<double> Tgrid(3);
	vector<double> rhogrid(2);
	Tgrid[0]   = 1.0;
	Tgrid[1]   = 2.0;
	Tgrid[2]   = 3.0;
	rhogrid[0] = 0.1;
	rhogrid[1] = 0.2;

	if ( match(gray->getTemperatureGrid(), Tgrid) )
	{
	    PASSMSG("Gray temperature grid correct.");
	}
	else
	{
	    FAILMSG("Gray temperature grid incorrect.");
	}

	if ( match(gray->getDensityGrid(), rhogrid) )
	{
	    PASSMSG("Gray density grid correct.");
	}
	else
	{
	    FAILMSG("Gray density grid incorrect.");
	}
    }

    // now reassign and check gray opacity for absorption opacities
    {
	gray = gray_abs;

	if (gray->getReactionType() != rtt_cdi::ABSORPTION) ITFAILS;
    }

    // Assign and check gray opacity 
    SP<rtt_cdi_test::DummyMultigroupOpacity> mg_total;
    SP<rtt_cdi_test::DummyMultigroupOpacity> mg_abs;
    mg_total = new rtt_cdi_test::DummyMultigroupOpacity();
    mg_abs   = new rtt_cdi_test::DummyMultigroupOpacity(rtt_cdi::ABSORPTION); 

    // check multigroup total opacities
    {
	mg = mg_total;

	if (!mg->data_in_tabular_form())             ITFAILS;
	if (mg->getReactionType() != rtt_cdi::TOTAL) ITFAILS;
	if (mg->getModelType() != rtt_cdi::ANALYTIC) ITFAILS;

	vector<double> Tgrid(3);
	vector<double> rhogrid(2);
	vector<double> egroups(4);
	Tgrid[0]   = 1.0;
	Tgrid[1]   = 2.0;
	Tgrid[2]   = 3.0;
	rhogrid[0] = 0.1;
	rhogrid[1] = 0.2;
	egroups[0] = 0.05;
	egroups[1] = 0.5;
	egroups[2] = 5.0;
	egroups[3] = 50.0;

	if ( match(mg->getTemperatureGrid(), Tgrid) )
	{
	    PASSMSG("Multigroup temperature grid correct.");
	}
	else
	{
	    FAILMSG("Multigroup temperature grid incorrect.");
	}

	if ( match(mg->getDensityGrid(), rhogrid) )
	{
	    PASSMSG("Multigroup density grid correct.");
	}
	else
	{
	    FAILMSG("Multigroup density grid incorrect.");
	}

	if ( match(mg->getGroupBoundaries(), egroups) )
	{
	    PASSMSG("Multigroup energy boundaries correct.");
	}
	else
	{
	    FAILMSG("Multigroup energy boundaries incorrect.");
	}

	if (mg->getNumTemperatures() != 3)    ITFAILS;
	if (mg->getNumDensities() != 2)       ITFAILS;
	if (mg->getNumGroupBoundaries() != 4) ITFAILS;
    }

    // noew reassign and check multigroup opacities for absorption
    {
	mg = mg_abs;

	if (mg->getReactionType() != rtt_cdi::ABSORPTION) ITFAILS;
    }
}

//---------------------------------------------------------------------------//

void gray_opacity_test()
{
    // ---------------------------- //
    // Create a GrayOpacity object. //
    // ---------------------------- //

    SP< GrayOpacity > spDGO;
    
    if ( spDGO = new rtt_cdi_test::DummyGrayOpacity() )
    {
	PASSMSG("SP to new GrayOpacity object created.");
    }
    else
    {
	FAILMSG("Unable to create a SP to new GrayOpacity object.");
    }

    // ------------------------ //
    // Dummy Gray Opacity Tests //
    // ------------------------ //
    
    double temperature          = 0.1;                          // keV
    double density              = 27.0;                         // g/cm^3
    double tabulatedGrayOpacity = temperature + density/1000.0; // cm^2/g
    
    double opacity = spDGO->getOpacity( temperature, density );

    if ( match ( opacity, tabulatedGrayOpacity ) ) 
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity computation was good.";
	PASSMSG( message.str() );
    }
    else
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity value is out of spec.";
	FAILMSG( message.str() );
    }

    // try using a vector of temps.

    std::vector< double > vtemperature(2);
    vtemperature[0] = 0.5;  // keV
    vtemperature[1] = 0.7;  // keV
    density         = 0.35; // g/cm^3
    std::vector< double > vRefOpacity( vtemperature.size() );
    for ( int i=0; i<vtemperature.size(); ++i )
	vRefOpacity[i] = vtemperature[i] + density/1000;

    std::vector< double > vOpacity = spDGO->getOpacity( vtemperature, 
							density );  

    if ( match ( vOpacity, vRefOpacity ) ) 
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity computation was good for a vector of temps.";
	PASSMSG( message.str() );
    }
    else
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity value is out of spec. for a vector of temps.";
	FAILMSG( message.str() );
    }
	

    // try using a vector of densities.

    std::vector< double > vdensity(5);
    vdensity[0] = 0.5;
    vdensity[1] = 1.0;
    vdensity[2] = 3.3;
    vdensity[3] = 5.0;
    vdensity[4] = 27.0;

    vRefOpacity.resize( vdensity.size() );
    for ( int i=0; i<vdensity.size(); ++i )
	vRefOpacity[i] = temperature + vdensity[i]/1000;

    vOpacity = spDGO->getOpacity( temperature, vdensity );

    if ( match ( vOpacity, vRefOpacity ) )
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity computation was good"
		<< " for a vector of densities.";
	PASSMSG( message.str() );
    }
    else
    {
	ostringstream message;
	message << spDGO->getDataDescriptor()
		<< " getOpacity value is out of spec."
		<< " for a vector of densities.";
	FAILMSG( message.str() );
    }
}

//---------------------------------------------------------------------------//

void multigroup_opacity_test()
{
    // ----------------------------------------- //
    // Create a Dummy Multigroup Opacity object. //
    // ----------------------------------------- //

    SP< MultigroupOpacity > spDmgO;
    
    if ( spDmgO = new rtt_cdi_test::DummyMultigroupOpacity() )
    {
	ostringstream message;
	message << "SP to new MultigroupOpacity object created.";
	PASSMSG( message.str() );
    }

    // --------------- //
    // MG Opacity test //
    // --------------- //
    
    // Setup the test point.
    double temperature = 0.01; // keV
    double density     = 2.0;  // g/cm^3

    // declare vector temps and densities
    vector<double> vtemperature(2);
    vector<double> vdensity(2);

    vtemperature[0] = 0.5; // keV
    vtemperature[1] = 0.7; // keV
    vdensity[0]     = 1.5; // g/cc
    vdensity[1]     = 2.0; // g/cc

    // The dummy opacity object should have 3 groups.  Check it.
    int ng = spDmgO->getNumGroupBoundaries()-1;
    if ( ng == 3 )
    {
	ostringstream message;
	message << "Correct number of groups found for "
		<< "MultigroupOpacity object.";
	PASSMSG( message.str() );
    }
    else
    {
	ostringstream message;
	message << "Wrong number of groups found for "
		<< "MultigroupOpacity object.";
	FAILMSG( message.str() );
    }

    const std::vector< double > energyBoundaries =
	spDmgO->getGroupBoundaries();

    // Create a container that hold all the MG opacities for a
    // specified temperature and density.  Fill this container with
    // the values that DummyMultigroupOpacity should contain.
    std::vector< double > tabulatedMGOpacity( ng );
    for ( int ig=0; ig<ng; ++ig)
	tabulatedMGOpacity[ig] = 2*(temperature + density/1000)
	    / (energyBoundaries[ig]+energyBoundaries[ig+1]);

    // Use the getOpacity accessor to obtain the MG opacities for a
    // specified temperature and density.
    std::vector< double > mgOpacity =
	spDmgO->getOpacity( temperature, density );

    // Make sure the accessor values match the expected values.
    if ( match( mgOpacity, tabulatedMGOpacity ) )
    {
	ostringstream message;
	message << spDmgO->getDataDescriptor()
		<< " getOpacity computation was good.";
	PASSMSG( message.str() );
    }
    else
    {
	ostringstream message;
	message << spDmgO->getDataDescriptor()
		<< " getOpacity value is out of spec.";
	FAILMSG( message.str() ); 
    }
    
    // Repeat with a vector of temps.
    
    // Reference values.

    // The opacity container is a vector<vector<double>>.  Each nested 
    // vector contains all of the group opacity values for a single
    // temperature. 

    // a MG opacity set for a single temperature, density combination
    // can be extracted from this container by using the following
    // type of assignment.
    // std::vector< double > vec1 = vRefMgOpacity[0];

    std::vector< std::vector< double > > vRefMgOpacity( ng );
    for ( int it=0; it<vtemperature.size(); ++it )
    {
	vRefMgOpacity[it].resize( ng );
	for ( int ig=0; ig<ng; ++ig )
	    vRefMgOpacity[it][ig] = 
		2.0 * ( vtemperature[it] + density/1000.0 ) 
		/ ( energyBoundaries[ig] + energyBoundaries[ig+1] );
    }
    
    // Retrieve the same set of opacity values via the getOpacity() accessor.
    std::vector< std::vector< double > > vMgOpacity 
	= spDmgO->getOpacity( vtemperature, density );
    
    // Compare the results.
    if ( match( vMgOpacity, vRefMgOpacity ) )
    {
	ostringstream message;
	message << spDmgO->getDataDescriptor()
		<< " getOpacity computation was good for a vector of  temps.";
	PASSMSG( message.str() );
    }
    else
    {
	ostringstream message;
	message << spDmgO->getDataDescriptor()
		<< " getOpacity value is out of spec. for a vector of temps.";
	FAILMSG( message.str() );
    }
    
    // STL-like accessor (MG opacities)

    // We have added STL-like getOpacity functions to DummyMultigroupOpacity,
    // these are not available through the rtt_cdi::MultigroupOpacity base
    // class so we test them as a DummyMultigroupOpacity.  This demonstrates
    // that one could make an opacity class that contains extra
    // functionality. Of course this functionality is not available through
    // CDI. 

    SP< rtt_cdi_test::DummyMultigroupOpacity > spDumMgOp;
    if ( spDumMgOp = new rtt_cdi_test::DummyMultigroupOpacity() )
    {
	ostringstream message;
	message << "SP to new DummyMultigroupOpacity object created.";
	PASSMSG( message.str() );
    }
    else
    {
	ostringstream message;
	message << "Unable to create a SP "
		<< "to a new DummyMultigroupOpacity object.";
	FAILMSG( message.str() );
    }

    // The STL-like accessors only work with 1-D containers.
    vector<double> vOpacity;
    vector<double> vRefOpacity;

    vOpacity.resize( vtemperature.size() * ng );
    vRefOpacity.resize( vtemperature.size() * ng );

    // Reference Values
    for ( int it=0; it<vtemperature.size(); ++it )
	for ( int ig=0; ig<ng; ++ig )
	    vRefOpacity[it*ng+ig] = 
		2.0 * ( vtemperature[it] + vdensity[it]/1000.0 )
		/ ( energyBoundaries[ig] + energyBoundaries[ig+1] );

    // Obtain values using getOpacity() accessor.
    spDumMgOp->getOpacity( vtemperature.begin(), vtemperature.end(),
			   vdensity.begin(), vdensity.end(),
			   vOpacity.begin() );
    
    // Compare the results:
    if ( match( vOpacity, vRefOpacity ) )
    {
	ostringstream message;
    	message << spDumMgOp->getDataDescriptor()
		<< " STL getOpacity() computation was good for a\n"
		<< " vector of temps. and a vector of densities.";
	PASSMSG( message.str() );
    }
    else
    {
	ostringstream message;
	message << spDumMgOp->getDataDescriptor()
		<< " STL getOpacity() value is out of spec. for a\n"
		<< " vector of temps. and a vector of densities.";
	FAILMSG( message.str() );
    }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_cdi::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	simple_tests();
	gray_opacity_test();
	multigroup_opacity_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tDummyOpacity, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_cdi_test::passed) 
    {
        cout << "**** tDummyOpacity Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tDummyOpacity." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tDummyOpacity.cc
//---------------------------------------------------------------------------//
