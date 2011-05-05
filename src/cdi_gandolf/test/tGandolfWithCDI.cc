//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/test/tGandolfWithCDI.cc
 * \author Thomas M. Evans
 * \date   Mon Oct 29 17:16:32 2001
 * \brief  Gandolf + CDI test.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "cdi_gandolf_test.hh"
#include "ds++/Release.hh"
#include "../GandolfFile.hh"
#include "../GandolfException.hh"
#include "../GandolfGrayOpacity.hh"
#include "../GandolfMultigroupOpacity.hh"
#include "cdi/CDI.hh" // this includes everything from CDI
#include "ds++/Assert.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

using rtt_cdi_gandolf::GandolfGrayOpacity;
using rtt_cdi_gandolf::GandolfMultigroupOpacity;
using rtt_cdi_gandolf::GandolfFile;
using rtt_cdi::GrayOpacity;
using rtt_cdi::MultigroupOpacity;
using rtt_cdi::CDI;
using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_gandolf_CDI()
{
	    
    // ----------------------------------------------- //
    // Test the data file "analyticOpacities.ipcress"  //
    // ----------------------------------------------- //
	    
    // -----------------------------------------------------------------
    // The Opacities in this file are computed from the following
    // analytic formula:
    //     opacity = rho * T^4,
    // rho is the density and T is the temperature.
    //
    // The grid in this data file has the following structure:
    //    T   = { 0.1, 1.0, 10.0 } keV.
    //    rho = { 0.1, 0.5, 1.0 } g/cm^3
    //    E_bounds = { 0.01, 0.03, 0.07, 0.1, 0.3, 0.7, 1.0, 3.0, 7.0 
    //                 10.0, 30.0, 70.0 100.0 } keV.
    //-----------------------------------------------------------------
	    
    // Gandolf data filename (IPCRESS format required)
    std::string op_data_file = "analyticOpacities.ipcress";
	    
    // ------------------------- //
    // Create GandolfFile object //
    // ------------------------- //
	    
    // Create a smart pointer to a GandolfFile object
    SP< const GandolfFile > spGFAnalytic;
	    
    // Try to instantiate the object.
    try 
    {
	spGFAnalytic = new const GandolfFile( op_data_file ); 
    }
    catch ( const rtt_cdi_gandolf::gmatidsException& GandError)
    {
	FAILMSG(GandError.what());
	ostringstream message;
	message << "Aborting tests because unable to instantiate "
		<< "GandolfFile object";
	FAILMSG(message.str());
	return;
    }
	    
    // If we make it here then spGFAnalytic was successfully instantiated.
    PASSMSG("SP to new GandolfFile object created (spGFAnalytic).");
			    
    // ----------------------------------- //
    // Create a GandolfGrayOpacity object. //
    // ----------------------------------- //
	    
    // Material identifier.  This data file has two materials: Al and
    // BeCu.  Al has the id tag "10001".
    const int matid=10001;
	    
    // Create a smart pointer to an opacity object.
    SP< const GrayOpacity > spOp_Analytic_ragray;
	    
    // Try to instantiate the opacity object.
    try
    {
	spOp_Analytic_ragray = new const GandolfGrayOpacity(
	    spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION ); 
    }
    catch ( const rtt_cdi_gandolf::GandolfException& GandError )
	// Alternatively, we could use:
	// catch ( rtt_cdi_gandolf::gkeysException GandError )
	// catch ( rtt_cdi_gandolf::gchgridsException GandError )
	// catch ( rtt_cdi_gandolf::ggetmgException GandError )
	// catch ( rtt_cdi_gandolf::ggetgrayException GandError )
    {
	ostringstream message;
	message << "Failed to create SP to new GandolfGrayOpacity object for "
		<< "Al_BeCu.ipcress data."
		<< std::endl << "\t" << GandError.what();
	FAILMSG(message.str());
	FAILMSG("Aborting tests.");
	return;
    }
	    
    // If we get here then the object was successfully instantiated.
    {
	ostringstream message;
	message << "SP to new GandolfGrayOpacity object created "
		<< "for analyticOpacities.ipcress.";	    
	PASSMSG(message.str());
    }

    // ----------------- //
    // Create CDI object //
    // ----------------- //
	    
    SP< CDI > spCDI_Analytic;
    if ( spCDI_Analytic = new CDI() )
    {
	ostringstream message;
	message << "SP to CDI object created successfully (GrayOpacity).";
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << "Failed to create SP to CDI object (GrayOpacity).";   
	FAILMSG(message.str()); 
    }    
	    
    // ------------------ //
    // Gray Opacity Tests //
    // ------------------ //

    // set the gray opacity
    spCDI_Analytic->setGrayOpacity(spOp_Analytic_ragray);
	    
    double temperature          = 10.0;                            // keV
    double density              = 1.0;                             // g/cm^3
    double tabulatedGrayOpacity = density * pow( temperature, 4 ); // cm^2/g

    rtt_cdi::Model    r = rtt_cdi::ROSSELAND;
    rtt_cdi::Reaction a = rtt_cdi::ABSORPTION;
	    
    double opacity = spCDI_Analytic->gray(r, a)->getOpacity(
	temperature, density );
	    
    if ( rtt_cdi_gandolf_test::match ( opacity, tabulatedGrayOpacity ) ) 
    {
	ostringstream message;
	message << spCDI_Analytic->gray(r, a)->getDataDescriptor()
		<< " getOpacity computation was good.";
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << spCDI_Analytic->gray(r, a)->getDataDescriptor()
		<< " getOpacity value is out of spec.";
	FAILMSG(message.str());
    }
	    
    // try using a vector of temps.
	    
    std::vector< double > vtemperature(2);
    vtemperature[0] = 0.5;  // keV
    vtemperature[1] = 0.7;  // keV
    density         = 0.35; // g/cm^3
    std::vector< double > vRefOpacity( vtemperature.size() );
    for ( size_t i=0; i<vtemperature.size(); ++i )
	vRefOpacity[i] = density * pow ( vtemperature[i], 4 );
	    
    std::vector< double > vOpacity = spCDI_Analytic->gray(r, a)->
	getOpacity( vtemperature, density );
	    
    if ( rtt_cdi_gandolf_test::match ( vOpacity, vRefOpacity ) ) 
    {
	ostringstream message;
	message << spCDI_Analytic->gray(r, a)->getDataDescriptor()
		<< " getOpacity computation was good for a vector of temps.";
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << spCDI_Analytic->gray(r, a)->getDataDescriptor()
		<< " getOpacity value is out of spec. for a vector of temps.";
	FAILMSG(message.str());
    }    

    // STL-like accessor
	    
    // The virtual base class does not support STL-like accessors
    // so we don't test this feature.
	    
    // Currently, KCC does not allow pure virtual + templates.
	    
	    
    // ----------------------------------------- //
    // Create a GandolfMultigorupOpacity object. //
    // ----------------------------------------- //
	    
    // Create a smart pointer to an opacity object.
    SP< const MultigroupOpacity > spOp_Analytic_ramg;
	    
    // Try to instantiate the opacity object.
    try
    {
	spOp_Analytic_ramg = new const GandolfMultigroupOpacity(
	    spGFAnalytic, matid, rtt_cdi::ROSSELAND, rtt_cdi::ABSORPTION);
    }
    catch ( const rtt_cdi_gandolf::GandolfException& GandError )
	// Alternatively, we could use:
	// catch ( rtt_cdi_gandolf::gkeysException GandError )
	// catch ( rtt_cdi_gandolf::gchgridsException GandError )
	// catch ( rtt_cdi_gandolf::ggetmgException GandError )
	// catch ( rtt_cdi_gandolf::ggetgrayException GandError )
    {
	ostringstream message;
	message << "Failed to create SP to new GandolfMultigroupOpacity "
		<< "object for Al_BeCu.ipcress data."
		<< std::endl << "\t" << GandError.what();
	FAILMSG(message.str());
	FAILMSG("Aborting tests.");
	return;
    }
	    
    // If we get here then the object was successfully instantiated.
    {
	ostringstream message;	
	message << "SP to new Gandolf multigroup opacity object created"
		<< "\n\tfor analyticOpacities.ipcress.";
	PASSMSG(message.str());
    }
	    
	    
    // ----------------------------------------------- //
    // Create a new CDI that has both Gray and MG data //
    // ----------------------------------------------- //

    // Add the multigroup opacity object to this CDI object.

    spCDI_Analytic->setMultigroupOpacity( spOp_Analytic_ramg );
	    
    // --------------- //
    // MG Opacity test //
    // --------------- //
	    
    // Set up the new test problem.
	    
    temperature = 0.3; // keV
    density     = 0.7; // g/cm^3
	    
    // This is the solution we compare against.
    int numGroups = 12;
    std::vector< double > tabulatedMGOpacity( numGroups );
    for ( int i=0; i<numGroups; ++i )
	tabulatedMGOpacity[i] = density * pow( temperature, 4 ); // cm^2/gm
	    
    // Request the multigroup opacity vector.
    std::vector< double > mgOpacity =
	spCDI_Analytic->mg(r, a)->getOpacity ( temperature, density );
	    
    if ( rtt_cdi_gandolf_test::match ( mgOpacity, tabulatedMGOpacity ) )
    {
	ostringstream message;
	message << spCDI_Analytic->mg(r, a)->getDataDescriptor()
		<< " getOpacity computation was good.";
	PASSMSG(message.str());
    }
    else
    {
	ostringstream message;
	message << spCDI_Analytic->mg(r, a)->getDataDescriptor()
		<< " getOpacity value is out of spec.";
	FAILMSG(message.str());
    }

    // Finally lets check to see if CDI catches some inappropriate accesses

    // multigroup access
    bool caught = false;
    try
    {
	spCDI_Analytic->mg(r, rtt_cdi::SCATTERING);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	PASSMSG("Good, caught illegal accessor to CDI-mg().");
	caught = true;
    }
    if (!caught)
	FAILMSG("Failed to catch illegal accessor to CDI-mg().");
    
    // gray access
    caught = false;
    try
    {
	spCDI_Analytic->gray(rtt_cdi::ANALYTIC, a);
    }
    catch (const rtt_dsxx::assertion &ass)
    {
	PASSMSG("Good, caught illegal accessor to CDI-gray().");
	caught = true;
    }
    if (!caught)
	FAILMSG("Failed to catch illegal accessor to CDI-gray().");
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
	test_gandolf_CDI();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tGandolfWithCDI, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_cdi_gandolf_test::passed) 
    {
        cout << "**** tGandolfWithCDI Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tGandolfWithCDI." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tGandolfWithCDI.cc
//---------------------------------------------------------------------------//
