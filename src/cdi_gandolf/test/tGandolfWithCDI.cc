//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tGandolfWithCDI.cc
 * \author Kelly Thompson
 * \date   Mon Mar 12 8:48:59 2001
 * \brief  Implementation file for tGandolfWithCDI
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>  // required by match()

#include "tGandolfWithCDI.hh"
#include "../Release.hh"

#include "../GandolfFile.hh"
#include "../GandolfException.hh"
#include "../GandolfGrayOpacity.hh"
#include "../GandolfGrayOpacity.t.hh"
#include "../GandolfMultigroupOpacity.hh"
#include "../GandolfMultigroupOpacity.t.hh"

#include "cdi/CDI.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

using rtt_cdi::ROSSELAND;
using rtt_cdi::ABSORPTION;

// --------------------- //
// Unit Test Frame Stuff //
// --------------------- //

namespace rtt_UnitTestFrame {
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
					   std::ostream& os_in ) {
	return rtt_dsxx::SP<TestApp> ( 
	    new rtt_gandolf_with_cdi_test::tGandolfWithCDI( argc, argv, os_in ));
    }
} // end namespace rtt_UnitTestFrame

// --------------------- //
// tGandolfWithCDI Stuff //
// --------------------- //

namespace rtt_gandolf_with_cdi_test 
{
    
    tGandolfWithCDI::tGandolfWithCDI( int argc, char *argv[], std::ostream& os_in )
	: rtt_UnitTestFrame::TestApp( argc, argv, os_in )
	{
	    os() << "Created tGandolfWithCDI." << std::endl;
	}
    
    std::string tGandolfWithCDI::version() const
	{
	    return rtt_cdi_gandolf::release();
	}
    
    
    //=======================================//
    // The tGandolfWithCDI tests start here. //
    //=======================================//
    
    
    std::string tGandolfWithCDI::runTest()
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
	    rtt_dsxx::SP< const rtt_cdi_gandolf::GandolfFile > spGFAnalytic;
	    
	    // Try to instantiate the object.
	    try 
		{
		    spGFAnalytic = new rtt_cdi_gandolf::GandolfFile( op_data_file ); 
		}
	    catch ( const rtt_cdi_gandolf::gmatidsException& GandError)
		{
		    fail() << std::endl << "\t" << GandError.what();
		    return "Unable to instantiate GandolfFile object spGFAnalytic.  Test sequence aborted.";
		}
	    
	    // If we make it here then spGFAnalytic was successfully instantiated.
	    pass() << "SP to new GandolfFile object created (spGFAnalytic).";
	    
	    
	    
	    // ----------------------------------- //
	    // Create a GandolfGrayOpacity object. //
	    // ----------------------------------- //
	    
	    // Material identifier.  This data file has two materials: Al and
	    // BeCu.  Al has the id tag "10001".
	    const int matid=10001;
	    
	    // Create a smart pointer to an opacity object.
	    rtt_dsxx::SP< const rtt_cdi::GrayOpacity > spOp_Analytic_ragray;
	    
	    // Try to instantiate the opacity object.
	    try
		{
		    spOp_Analytic_ragray 
			= new rtt_cdi_gandolf::GandolfGrayOpacity(
			    spGFAnalytic,
			    matid,
			    rtt_cdi::ROSSELAND,         // enumeration
			    rtt_cdi::ABSORPTION );      // enumeration
		}
	    catch ( const rtt_cdi_gandolf::GandolfException& GandError )
		// Alternatively, we could use:
		// catch ( rtt_cdi_gandolf::gkeysException GandError )
		// catch ( rtt_cdi_gandolf::gchgridsException GandError )
		// catch ( rtt_cdi_gandolf::ggetmgException GandError )
		// catch ( rtt_cdi_gandolf::ggetgrayException GandError )
		{
		    fail() << "Failed to create SP to new GandolfGrayOpacity object for "
			   << "Al_BeCu.ipcress data."
			   << std::endl << "\t" << GandError.what();
		    return "Unable to instantiate GandolfGrayOpacity object.  Test sequence aborted.";
		}
	    
	    // If we get here then the object was successfully instantiated.
	    pass() << "SP to new GandolfGrayOpacity object created for analyticOpacities.ipcress.";
	    
	    // ----------------- //
	    // Create CDI object //
	    // ----------------- //
	    
	    rtt_dsxx::SP< rtt_cdi::CDI > spCDI_Analytic;
	    if ( spCDI_Analytic = new rtt_cdi::CDI() )
		pass() << "SP to CDI object created successfully (GrayOpacity).";
	    else
		fail() << "Failed to create SP to CDI object (GrayOpacity).";
	    
	    
	    // ------------------ //
	    // Gray Opacity Tests //
	    // ------------------ //

	    // set the gray opacity
	    spCDI_Analytic->setGrayOpacity(spOp_Analytic_ragray);
	    
	    double temperature = 10.0; // keV
	    double density = 1.0; // g/cm^3
	    double tabulatedGrayOpacity = density * pow( temperature, 4 ); // cm^2/g
	    
	    double opacity = spCDI_Analytic->gray(ROSSELAND,ABSORPTION)->getOpacity( temperature, density );
	    
	    if ( match ( opacity, tabulatedGrayOpacity ) ) 
		pass() << spCDI_Analytic->gray(ROSSELAND,ABSORPTION)->getDataDescriptor()
		       << " getOpacity computation was good.";
	    else
		fail() << spCDI_Analytic->gray(ROSSELAND,ABSORPTION)->getDataDescriptor()
		       << " getOpacity value is out of spec.";
	    
	    // try using a vector of temps.
	    
	    std::vector< double > vtemperature(2);
	    vtemperature[0] = 0.5; // keV
	    vtemperature[1] = 0.7; // keV
	    density = 0.35; // g/cm^3
	    std::vector< double > vRefOpacity( vtemperature.size() );
	    for ( int i=0; i<vtemperature.size(); ++i )
		vRefOpacity[i] = density * pow ( vtemperature[i], 4 );
	    
	    std::vector< double > vOpacity = spCDI_Analytic->gray(ROSSELAND,ABSORPTION)->
		getOpacity( vtemperature, density );
	    
	    if ( match ( vOpacity, vRefOpacity ) ) 
		pass() << spCDI_Analytic->gray(ROSSELAND,ABSORPTION)->getDataDescriptor()
		       << " getOpacity computation was good for a vector of temps.";
	    else
		fail() << spCDI_Analytic->gray(ROSSELAND,ABSORPTION)->getDataDescriptor()
		       << " getOpacity value is out of spec. for a vector of temps.";
	    
	    
	    // STL-like accessor
	    
	    // The virtual base class does not support STL-like accessors
	    // so we don't test this feature.
	    
	    // Currently, KCC does not allow pure virtual + templates.
	    
	    
	    // ----------------------------------------- //
	    // Create a GandolfMultigorupOpacity object. //
	    // ----------------------------------------- //
	    
	    // Create a smart pointer to an opacity object.
	    rtt_dsxx::SP< const rtt_cdi::MultigroupOpacity > spOp_Analytic_ramg;
	    
	    // Try to instantiate the opacity object.
	    try
		{
		    spOp_Analytic_ramg
			= new rtt_cdi_gandolf::GandolfMultigroupOpacity(
			    spGFAnalytic,
			    matid,
			    rtt_cdi::ROSSELAND,         // enumeration
			    rtt_cdi::ABSORPTION );      // enumeration
		}
	    catch ( const rtt_cdi_gandolf::GandolfException& GandError )
		// Alternatively, we could use:
		// catch ( rtt_cdi_gandolf::gkeysException GandError )
		// catch ( rtt_cdi_gandolf::gchgridsException GandError )
		// catch ( rtt_cdi_gandolf::ggetmgException GandError )
		// catch ( rtt_cdi_gandolf::ggetgrayException GandError )
		{
		    fail() << "Failed to create SP to new GandolfMultigroupOpacity "
			   << "object for Al_BeCu.ipcress data."
			   << std::endl << "\t" << GandError.what();
		    return "Unable to instantiate GandolfMultigroupOpacity object.  Test sequence aborted.";
		}
	    
	    // If we get here then the object was successfully instantiated.
	    pass() << "SP to new Gandolf multigroup opacity object created"
		   << "\n\tfor analyticOpacities.ipcress.";
	    
	    
	    // ----------------------------------------------- //
	    // Create a new CDI that has both Gray and MG data //
	    // ----------------------------------------------- //
	    
// 	    if ( spCDI_Analytic = new rtt_cdi::CDI( spOp_Analytic_ragray ) ) 
// 		pass() << "SP to CDI object created successfully (GrayOpacity + Multigroup).";
// 	    else
// 		fail() << "Failed to create SP to CDI object (GrayOpacity + Multigroup).";

	    // Add a multigroup opacity object to this CDI object.

	    spCDI_Analytic->setMultigroupOpacity( spOp_Analytic_ramg );
	    
	    // --------------- //
	    // MG Opacity test //
	    // --------------- //
	    
	    // Set up the new test problem.
	    
	    temperature = 0.3; // keV
	    density = 0.7; // g/cm^3
	    
	    // This is the solution we compare against.
	    int numGroups = 12;
	    std::vector< double > tabulatedMGOpacity( numGroups );
	    for ( int i=0; i<numGroups; ++i )
		tabulatedMGOpacity[i] = density * pow( temperature, 4 ); // cm^2/gm
	    
	    // Request the multigroup opacity vector.
	    std::vector< double > mgOpacity =
		spCDI_Analytic->mg(ROSSELAND,ABSORPTION)->getOpacity ( temperature, density );
	    
	    if ( match ( mgOpacity, tabulatedMGOpacity ) )
		pass() << spCDI_Analytic->mg(ROSSELAND,ABSORPTION)->getDataDescriptor()
		       << " getOpacity computation was good.";
	    else
		fail() << spCDI_Analytic->mg(ROSSELAND,ABSORPTION)->getDataDescriptor()
		       << " getOpacity value is out of spec.";
	    
	    
	    // ---------------------- //
	    // Print the test result. //
	    // ---------------------- //
	    
	    if ( passed() ) {
		pass() << "All tests passed.";
		return "All tests passed.";
	    }
	    
	    return "Some tests failed.";
	    
	} // end of runTest()
    
    
    
    //---------------------------------------------
    // Compare Reference value to computed values
    //---------------------------------------------
    
    bool tGandolfWithCDI::match( const double computedValue,
				 const double referenceValue ) const
	{
	    // Compare items up to 10 digits of accuracy.
	    const double TOL = 1.0e-10;
	    
	    // Calculate the absolute value of the relative difference between 
	    // the computed and reference values.
	    double reldiff = fabs( ( computedValue - referenceValue )
				   / referenceValue );
	    
	    // If the comparison fails then return "false" to indicate that
	    // the test failed.
	    if ( reldiff > TOL ) return false;
	    
	    return true;    
	    
	} // end of tGandolfWithCDI::match( double, double )
    
    bool tGandolfWithCDI::match( 
	const std::vector< double >& computedValue, 
	const std::vector< double >& referenceValue ) const
	{
	    // If the vector sizes don't match then die
	    if ( computedValue.size() != referenceValue.size() )
		return false;
	    
	    // Compare items up to 10 digits of accuracy.
	    const double TOL = 1.0e-10;
	    
	    // Test each item in the list
	    double reldiff = 0.0;
	    for ( int i=0; i<computedValue.size(); ++i )
		{
		    reldiff = fabs( ( computedValue[i] - referenceValue[i] )
				    / referenceValue[i] );
		    // If the comparison fails then return "false" to indicate 
		    // that the test failed.
		    if ( reldiff > TOL ) return false;
		}
	    return true;
	} // end of tGandolfWithCDI::match( vector<double>, vector<double> )
    
    bool tGandolfWithCDI::match( 
	const std::vector< std::vector< double > >& computedValue, 
	const std::vector< std::vector< double > >& referenceValue ) const 
	{
	    // Compare items up to 10 digits of accuracy.
	    const double TOL = 1.0e-10;
	    
	    // Test each item in the list
	    double reldiff = 0.0;
	    for ( int i=0; i<computedValue.size(); ++i )
		{
		    for ( int j=0; j<computedValue[i].size(); ++j )
			{
			    reldiff = fabs( ( computedValue[i][j] - referenceValue[i][j] )
					    / referenceValue[i][j] );
			    // If the comparison fails then stop testing and
			    // return "false" to indicate that the test
			    // failed. 
			    if ( reldiff > TOL ) return false; 
			}
		}
	    return true;
	} // end of tGandolfWithCDI::match( vector<double>, vector<double> )
    
} // end namespace rtt_gandolf_with_cdi_test

//---------------------------------------------------------------------------//
// end of tGandolfWithCDI.cc
//---------------------------------------------------------------------------//
