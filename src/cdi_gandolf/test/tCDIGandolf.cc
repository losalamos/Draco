//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tCDIGandolf.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 13:07:00 2000
 * \brief  Implementation file for tCDIGandolf
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "tCDIGandolf.hh"
#include "../GandolfFile.hh"
#include "../GandolfOpacity.hh"
#include "../GandolfException.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"
//#include "cdi/Opacity.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// Unit Test Frame Stuff
//----------------------------------------
namespace rtt_UnitTestFrame {
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
					   std::ostream& os_in ) {
	return rtt_dsxx::SP<TestApp> ( 
	    new rtt_cdi_gandolf_test::tCDIGandolf( argc, argv, os_in ));
    }
} // end namespace rtt_UnitTestFrame

// tCDIGandolf Stuff
//--------------------------------------------------
namespace rtt_cdi_gandolf_test {

using std::string;
using std::cout;
using std::endl;
using std::vector;
using rtt_dsxx::SP;

tCDIGandolf::tCDIGandolf( int argc, char *argv[], std::ostream& os_in )
    : rtt_UnitTestFrame::TestApp( argc, argv, os_in )
{
    os() << "Created tCDIGandolf" << endl;
}

string tCDIGandolf::version() const
{
    return rtt_cdi_gandolf::release();
}

string tCDIGandolf::runTest()
{
    // Gandolf data filename (IPCRESS format required)
    string op_data_file = "Al_BeCu.ipcress";
    
    // Material identifier.  This data file has two materials: Al and
    // BeCu.  Al has the id tag "10001".
    const int matid=10001;
    
    // ------------------------- //
    // Create GandolfFile object //
    // ------------------------- //
    
    SP<rtt_cdi_gandolf::GandolfFile> spGFABC;

    // Attempt to instantiate the object.
    try
	{
	    spGFABC = new rtt_cdi_gandolf::GandolfFile( op_data_file ); 
	}
    catch ( rtt_cdi_gandolf::GandolfException gerr )
	{
	    fail() << std::endl << "\t" << gerr.errorSummary();
	    return "Some tests failed.";
	}

    // If we make it here then spGFABC was successfully instantiated.
    pass() << "SP to new GandolfFile object created for Al_BeCu.ipcress data.";

    // Test the GandolfFile object.
    if ( spGFABC->getDataFilename() == op_data_file )
	pass() << "GandolfFile object is now linked to the Al_BeCu.ipcress data file.";
    else
	fail() << "GandolfFile object failed to link itself to theAl_BeCu.ipcress  data file.";

    if ( spGFABC->getNumMaterials() == 2 )
	pass() << "The correct number of materials was found in the Al_BeCu.ipcress data file.";
    else
	fail() << "spGFABC did not find the correct number of materials in the Al_BeCu.ipcress data file.";

    // ---------------------- //
    // Create Opacity object. //
    // ---------------------- //

    // Create a smart pointer to an Opacity object.
    SP<rtt_cdi::Opacity> spOpacityABC;

    // Try to instantiate the Opacity object.
    try 
	{
	    spOpacityABC 
		= new rtt_cdi_gandolf::GandolfOpacity( spGFABC, matid );
	}
    catch ( rtt_cdi_gandolf::GandolfException gerr )
	// Alternatively, we could use:
	// catch ( rtt_cdi_gandolf::gkeysException gerr )
	{
	    fail() << "Failed to create SP to new Opacity object for Al_BeCu.ipcress data."
		   << std::endl << "\t" << gerr.errorSummary();
	    return "Some tests failed.";
	}

    // If we get here then the object was successfully instantiated.
    pass() << "SP to new Opacity object created for Al_BeCu.ipcress data.";

    // ----------------- //
    // Gray Opacity Test //
    // ----------------- //

    double temperature = 0.1; // keV
    double density = 27.0; // g/cm^3
    double tabulatedGrayOpacity = 4271.7041147070677; // cm^2/g

    if ( ! testGrayRosselandOpacityAccessorPassed( 
	spOpacityABC, temperature, density, tabulatedGrayOpacity ) )
	return "testGrayRosselandOpacityAccessorPassed() failed.";
    
  // --------------- //
  // MG Opacity test //
  // --------------- //
  
    temperature = 0.01; // keV
    density = 2.0; // g/cm^3

    // The solution to compare against:
    int numGroups = 33;
    double tabulatedMGOpacityArray[] =
    {   2.4935160506e+08, 
	2.6666668566e+04, 
	1.6270311515e+04, 
	1.7634614752e+04, 
	4.4986543915e+04, 
	9.9917288352e+04, 
	8.3261040205e+04, 
	5.9742699269e+04, 
	4.0372986257e+04, 
	2.6156310247e+04, 
	1.6356517677e+04, 
	1.0006990734e+04, 
	5.9761421730e+03, 
	3.5201139934e+03, 
	2.0761987266e+03, 
	6.8547729207e+03, 
	4.1253621873e+03, 
	2.4194579352e+03, 
	1.3888963339e+03, 
	7.8973170438e+02, 
	4.3992670770e+02, 
	2.4390304515e+02, 
	1.3380845641e+02, 
	6.9706128602e+01,
	3.7034400058e+01, 
	1.9776964848e+01, 
	1.0437495354e+01, 
	5.4952246461e+00, 
	2.8823716869e+00, 
	1.4965372557e+00, 
	7.8893806873e-01, 
	4.1875736190e-01, 
	2.1668498635e-01 }; // KeV, numGroups entries.
     vector<double> tabulatedMGOpacity(numGroups);
     std::copy( tabulatedMGOpacityArray, 
		tabulatedMGOpacityArray+numGroups,
		tabulatedMGOpacity.begin() );    

     if ( ! testMGRosselandOpacityAccessorPassed( 
	 spOpacityABC, temperature, density, tabulatedMGOpacity ) )
	 return "testGrayRosselandOpacityAccessorPassed() failed.";
     
     // Interpolate the multigroup opacities for T = 0.01 keV and
     // density = 2.0 g/cm^3. 
     vector<double> mgRosselandOpacity;
     try
	 {
	     // The tabulated data is for absoption only so we over
	     // ride the default total opacity with rosseland
	     // absorption opacities
	     mgRosselandOpacity	
		 = spOpacityABC->getMGRosseland( 0.01, 2.0, "ramg" );
	 }
     catch ( rtt_cdi_gandolf::ggetmgException gerr )
	 {
	    fail() << std::endl << "\t" << gerr.errorSummary();
	    return "Some tests failed.";
	 }

     // Compare the interpolated value with previous interpolations:
     if ( match( mgRosselandOpacity, tabulatedMGOpacity ) )
 	pass() << "Multigroup Opacity computation was good for Al_BeCu.ipcress data.";
     else
 	fail() << "Multigroup Opacity computation failed for Al_BeCu.ipcress data.";

// ----------------------------------------------- //
// Test the data file "analyticOpacities.ipcress"  //
// ----------------------------------------------- //

     // The Opacities in this file are computed from the following
     // analytic formula:
     //     opacity = rho * T^4,
     // rho is the density and T is the temperature.

     // Gandolf data filename (IPCRESS format required)
     op_data_file = "analyticOpacities.ipcress";

     // ------------------------- //
     // Create GandolfFile object //
     // ------------------------- //

     // Create a smart pointer to a GandolfFile object
     SP<rtt_cdi_gandolf::GandolfFile> spGFAnalytic;

     // Try to instantiate the object.
     try 
	 {
	     spGFAnalytic = new rtt_cdi_gandolf::GandolfFile( op_data_file ); 
	 }
     catch ( rtt_cdi_gandolf::gmatidsException gerr)
	 {
	     fail() << std::endl << "\t" << gerr.errorSummary();
	     return "Some tests failed.";
	 }
     
     // If we make it here then spGFAnalytic was successfully instantiated.
     pass() << "SP to new GandolfFile object created (spGFAnalytic).";
     
     // Test the GandolfFile object.
     if ( spGFAnalytic->getDataFilename() == op_data_file )
	 pass() << "GandolfFile object is now linked to the data file.";
     else
	 fail() << "GandolfFile object failed to link itself to the data file.";

     if ( spGFAnalytic->getNumMaterials() == 1 )
	 pass() << "The correct number of materials was found in the data file.";
     else
	 fail() << "spGFAnalytic did not find the correct number of materials in the data file.";

     // --------------------- //
     // Create Opacity object //
     // --------------------- //

     // Create a smart pointer to an Opacity object.
     SP<rtt_cdi::Opacity> spOpacityAnalytic;
     
     // Try to instantiate the Opacity object.
     try 
	 {
	     spOpacityAnalytic 
		 = new rtt_cdi_gandolf::GandolfOpacity( spGFAnalytic, matid );
	 }
     catch ( rtt_cdi_gandolf::gkeysException gerr )
	 {
	    fail() << std::endl << "\t" << gerr.errorSummary();
	    return "Some tests failed.";
	 }
     catch ( rtt_cdi_gandolf::gchgridsException gerr )
	 {
	    fail() << std::endl << "\t" << gerr.errorSummary();
	    return "Some tests failed.";
	 }

     // If we get here then the object was successfully instantiated.
     pass() << "SP to new Opacity object created for analyticOpacities.ipcress.";

     // ----------------- //
     // Gray Opacity Test //
     // ----------------- //

     temperature = 1.0; // keV
     density = 1.0; // g/cm^3
     
     // The analytic formula shown above does not include the
     // addition of scattering. 
     double scatter = 0.4; // cm^2/g
     tabulatedGrayOpacity = density * pow( temperature, 4 ) + scatter;

     if ( ! testGrayRosselandOpacityAccessorPassed( 
	 spOpacityAnalytic, temperature, density, tabulatedGrayOpacity ) )
	 return "testGrayRosselandOpacityAccessorPassed() failed.";
     
     //---------------- //
     // MG Opacity test //
     //---------------- //
     
     temperature = 0.3; // keV
     density = 0.7; // g/cm^3
     
     // This is the solution we compare against.
     numGroups = 12;
     tabulatedMGOpacity.resize(12);
     for ( int i=0; i<numGroups; ++i )
	 tabulatedMGOpacity[i] = density * pow( temperature, 4 ); // cm^2/gm

     if ( ! testMGRosselandOpacityAccessorPassed( 
	 spOpacityAnalytic, temperature, density, tabulatedMGOpacity ) )
	 return "testMGRosselandOpacityAccessorPassed() failed.";
     
    // ------------------------------------------------------------ //
    // Test the Plank routines using analyticOpacities.ipcress data //
    // ------------------------------------------------------------ //
    
    // The Opacities in this file are computed from the following
    // analytic formula:
    //     opacity = rho * T^4,
    // rho is the density and T is the temperature.
    
    // spGFAnalytic already points to the correct file so we don't repeat the 
    // coding.
    
    // Dito for spOpacityAnalytic.
    
    // ----------------- //
    // Gray Opacity Test //
    // ----------------- //
    
    if ( ! testGrayPlankOpacityAccessorPassed( spOpacityAnalytic ) )
	return "testGrayPlankOpacityAccessor() failed.";
    
    // --------------- //
    // MG Opacity test //
    // --------------- //
    
    // If this test fails then stop testing.
    if ( ! testMGPlankOpacityAccessorPassed( spOpacityAnalytic ) )
	return "testMGPlankOpacityAccessor() failed.";
    
    // ------------------------ //
    // Access temperature grid. //
    // ------------------------ //
    
    testTemperatureGridAccessor( spOpacityAnalytic );
    
    // ------------------------ //
    // Access the density grid. //
    // ------------------------ //
    
    testDensityGridAccessor( spOpacityAnalytic );
    
    // ----------------------------- //
    // Access the energy boundaries. //
    // ----------------------------- //
    
    testEnergyBoundaryAccessor( spOpacityAnalytic );
    
    // ----------------------
    // Print the test result.
    // ----------------------
    if (passed()) {
	pass() << "All tests passed.";
	return "All tests passed.";
    }
    return "Some tests failed.";
}

//---------------------------------------------
// Compare Reference value to computed values
//---------------------------------------------

bool tCDIGandolf::match( vector<double> computedValue, 
			 vector<double> referenceValue )
{
    // Start by assuming that the two quantities match exactly.
    bool em = true;

    // Compare items up to 10 digits of accuracy.
    const double TOL = 1.0e-10;

    // Test each item in the list
    double absdiff = 0.0;
    for ( int i=0; i<computedValue.size(); ++i )
	{
	    absdiff = fabs( ( computedValue[i] - referenceValue[i] )
			    / referenceValue[i] );
	    // If the comparison fails then change the value of "em"
	    // and exit the loop.
	    if ( absdiff > TOL )
		{
		    em = false;
		    break;
		}
	}
    return em;
} // end of tCDIGandolf::match( vector<double>, vector<double> )

bool tCDIGandolf::match( const double computedValue,
			 const double referenceValue )
{
    // Start by assuming that the two quantities match exactly.
    bool em = true;

    // Compare items up to 10 digits of accuracy.
    const double TOL = 1.0e-10;

    // Calculate the absolute value of the relative difference between 
    // the computed and reference values.
    double absdiff = fabs( ( computedValue - referenceValue )
			   / referenceValue );
    
    // If the comparison fails then change the value of "em" return
    // the result;
    if ( absdiff > TOL )
	em = false;

    return em;    

} // end of tCDIGandolf::match( double, double )

// ---------------------------------------- //
// Test the gray Rosseland opacity accessor // 
// ---------------------------------------- //

bool tCDIGandolf::testGrayRosselandOpacityAccessorPassed(
    const SP<rtt_cdi::Opacity> spOpacity, 
    const double temperature, const double density, 
    const double tabulatedValue )
{
    double grayRosselandOpacity;
    try 
	{
	    grayRosselandOpacity 
		= spOpacity->getGrayRosseland( temperature, density );
	}
    catch ( rtt_cdi_gandolf::ggetgrayException gerr )
	{
	    fail() << std::endl << "\t" << gerr.errorSummary();
	    return false;
	    //return "Some tests failed.";
	}

    // Make sure that the interpolated value matches previous
    // interpolations. 

    if ( match( grayRosselandOpacity, tabulatedValue ) )
	pass() << "grayRosselandOpacity computation was good for "
	       << spOpacity->getDataFilename(); 
    else
	fail() << "grayRosselandOpacity value is out of spec. for "
	       << spOpacity->getDataFilename();
    
    return true;
}

// ---------------------------------------------- //
// Test the multigroup Rosseland opacity accessor //
// ---------------------------------------------- //

bool tCDIGandolf::testMGRosselandOpacityAccessorPassed(
    const SP<rtt_cdi::Opacity> spOpacity,
    const double temperature, const double density,
    const vector<double> tabulatedValues )
{
    // Interpolate the multigroup opacities for T = 0.3 keV and
    // density = 0.7 g/cm^3. 
    vector<double> mgRosselandOpacity( tabulatedValues.size() );
    try
	{
	    // override default "rtmg" mode for analytic values.
	    mgRosselandOpacity
		= spOpacity->getMGRosseland( temperature, density, "ramg" );
	}
    catch ( rtt_cdi_gandolf::ggetmgException gerr )
	{
	    fail() << std::endl << "\t" << gerr.errorSummary();
	    return false;
	    //return "Some tests failed.";
	}

     // print the results
//      cout << endl << "Multigroup Rosseland Opacities for Aluminum (matID=" 
// 	  << matid << ") test:" << endl
// 	  << "   MGOpacities at T=0.3 keV and rho=0.7 g/cm^3 ( in cm^2/gm. )."
// 	  << endl << endl;
//      for ( int i=0; i<mgRosselandOpacity.size(); ++i )
// 	 cout << "   Opacity(group=" << std::setw(2) << i << ") = " 
// 	      << std::scientific << std::setprecision(5) 
// 	      << mgRosselandOpacity[i] << endl;
//      cout << endl;
    

//     // Compare the interpolated value with previous interpolations:
     if ( match( mgRosselandOpacity, tabulatedValues ) )
 	pass() << "Multigroup Rosseland Opacity computation was\n\t"
	       << "good for the data obtained from the file " 
	       << "\"" << spOpacity->getDataFilename() << "\".";
     else
 	fail() << "Multigroup Rosseland Opacity computation failed for \n\t"
	       << "the data obtained from the file" 
	       << "\"" << spOpacity->getDataFilename() << "\".";

     return true;

} // end of tCDIGandolf::testGrayRosselandOpacityAccessorPassed( )

// ------------------------------------ //
// Test the gray Plank opacity accessor //
// ------------------------------------ //

bool tCDIGandolf::testGrayPlankOpacityAccessorPassed(
    SP<rtt_cdi::Opacity> spOpacity )
{
    // Obtain a Plank Gray Opacity value for T=3.0 keV and density 
    // = 0.7 g/cm^3.
    double T = 3.0;   // keV
    double rho = 0.7; // g/cm^3
    double grayPlankOpacity = 0.0;
    try 
	{
	    grayPlankOpacity 
		= spOpacity->getGrayPlank( T, rho );
	}
    catch ( rtt_cdi_gandolf::ggetgrayException gerr )
	{
	    fail() << std::endl << "\t" << gerr.errorSummary();
	    return false;
	    // return "Some tests failed.";
	}
     
    // Make sure that the interpolated value matches previous
    // interpolations. 
    double tabulatedGrayOpacity = rho*pow(T,4); // cm^2/g
    
    if ( match( grayPlankOpacity, tabulatedGrayOpacity ) )
	pass() << "grayPlankOpacity computation was good for analyticOpacity data.";
    else
	fail() << "grayPlankOpacity value is out of spec for analyticOpacity data.";
    
    return true;

} // end of bool tCDIGandolf::testGrayOpacityAccessorPassed( )

// ------------------------------------------ //
// Test the multigroup Plank Opacity accessor //
// ------------------------------------------ //

bool tCDIGandolf::testMGPlankOpacityAccessorPassed( 
    SP<rtt_cdi::Opacity> spOpacity )
{
    // The solution to compare against:
    int numGroups = 12;
    vector<double> tabulatedMGOpacity(numGroups);
    
    double T   = 0.4; // keV
    double rho = 0.22; // g/cm^3
    for ( int i=0; i<numGroups; ++i )
	tabulatedMGOpacity[i] = rho*pow(T,4); // cm^2/gm
    
    vector<double> mgPlankOpacity;
    
    // Interpolate the multigroup opacities for T = 0.4 keV and
    // density = 2.2 g/cm^3. 
    try
	{
	    mgPlankOpacity
		= spOpacity->getMGPlank( T, rho );
	}
    catch ( rtt_cdi_gandolf::ggetmgException gerr )
	{
	    fail() << std::endl << "\t" << gerr.errorSummary();
	    // return "Some tests failed.";
	    return false;
	}
    
     // print the results
//      cout << endl << "Multigroup Plank Opacities for Aluminum (matID=" 
// 	  << matid << ") test:" << endl
// 	  << "   MGOpacities at T=0.4 keV and rho=0.22 g/cm^3 ( in cm^2/gm. )."
// 	  << endl << endl;
//      for ( int i=0; i<mgPlankOpacity.size(); ++i )
// 	 cout << "   Opacity(group=" << std::setw(2) << i << ") = " 
// 	      << std::scientific << std::setprecision(5) 
// 	      << mgPlankOpacity[i] << endl;
//      cout << endl;

    // Compare the interpolated value with previous interpolations:
    if ( match( mgPlankOpacity, tabulatedMGOpacity ) )
	pass() << "Multigroup Plank Opacity computation was good "
	       << "for analyticOpacity data.";
    else
	fail() << "Multigroup Plank Opacity computation failed "
	       << "for analyticOpacity data.";

    return true;

} // end of void tCDIGandolf::testMGOpacityAccessor( )

// ---------------------------------- //
// Test the Temperature Grid Accessor //
// ---------------------------------- //

void tCDIGandolf::testTemperatureGridAccessor( 
    SP<rtt_cdi::Opacity> spOpacity )
{
     // Read the temperature grid from the IPCRESS file.     
     vector<double> temps = spOpacity->getTemperatureGrid();

     // Verify that the size of the temperature grid looks right.  If
     // it is the right size then compare the temperature grid data to 
     // the data specified when we created the IPCRESS file using TOPS.
     if ( temps.size() == spOpacity->getNumTemperatures() &&
	  temps.size() == 3 )
	 {
	     pass() << "The number of temperature points found in the data\n\t" 
		    << "grid matches the number returned by the\n\t"
		    << "getNumTemperatures() accessor.";

	     // The grid specified by TOPS has 3 temperature points.
	     vector<double> temps_ref( temps.size() );
	     temps_ref[0] = 0.1;
	     temps_ref[1] = 1.0;
	     temps_ref[2] = 10.0;

	     // Compare the grids.
	     if ( match( temps, temps_ref ) )
		 pass() << "Temperature grid matches.";
	     else
		 fail() << "Temperature grid did not match.";
	     
	     // cout << "Temperature Grid (keV):" << endl << endl;
	     // for ( int i=0; i<temps.size(); ++i )
	     //   cout << "   T[" << i << "] = " << temps[i] << endl;
	     
	 }
     else
	 {
	     fail() << "The number of temperature points found in the data\n\t"
		    << "grid does not match the number returned by the\n\t"
		    << "getNumTemperatures() accessor.";
	     fail() << "Did not test the results returned by\n\t"
		    << "getTemperatureGrid().";
	 }

} // end of  tCDIGandolf::testTemperatureGridAccessor( )

// ------------------------------ //
// Test the Density Grid Accessor //
// ------------------------------ //

void tCDIGandolf::testDensityGridAccessor( 
    SP<rtt_cdi::Opacity> spOpacity )
{
    // Read the grid from the IPCRESS file.     
    vector<double> density = spOpacity->getDensityGrid();
    
    // Verify that the size of the density grid looks right.  If
    // it is the right size then compare the density grid data to 
    // the data specified when we created the IPCRESS file using TOPS.
    if ( density.size() == 3 &&
	 density.size() == spOpacity->getNumDensities() )
	{
	    pass() << "The number of density points found in the data\n\t"
		   << "grid matches the number returned by the\n\t"
		   << "getNumDensities() accessor.";
	    
	    // The grid specified by TOPS has 3 density points
	    vector<double> density_ref( density.size() );
	    density_ref[0] = 0.1;
	    density_ref[1] = 0.5;
	    density_ref[2] = 1.0;
	    
	    // Compare the grids.
	    if ( match( density, density_ref ) )
		pass() << "Density grid matches.";
	    else
		fail() << "Density grid did not match.";
	    
// 	     cout << "Density Grid (g/cm^3):" << endl << endl;
// 	     for ( int i=0; i<density.size(); ++i )
// 		 cout << "   rho[" << i << "] = " << density[i] << endl;
	}
    else
	{
	    fail() << "The number of density points found in the data\n\t"
		   << "grid does not match the number returned by the\n\t"
		   << "getNumDensities() accessor.";
	    fail() << "Did not test the results returned by\n\t"  
		   << "getDensityGrid().";
	}
} // end of void tCDIGandolf::testDensityGridAccessor( )

// --------------------------------- //
// Test the Energy Boundary Accessor //
// --------------------------------- //

void tCDIGandolf::testEnergyBoundaryAccessor( SP<rtt_cdi::Opacity> spOpacity )
{
    
     // Read the grid from the IPCRESS file.     
     vector<double> ebounds = spOpacity->getGroupBoundaries();

     // Verify that the size of the group boundary grid looks right.  If
     // it is the right size then compare the energy groups grid data to 
     // the data specified when we created the IPCRESS file using TOPS.
     if ( ebounds.size() == 13 &&
	  ebounds.size() == spOpacity->getNumGroupBoundaries() )
	 {
	     pass() << "The number of energy boundary points found in the data\n\t"
		    << "grid matches the number returned by the\n\t"
		    << "getNumGroupBoundaries() accessor.";

	     // The grid specified by TOPS has 13 energy boundaries.
	     vector<double> ebounds_ref(ebounds.size());
	     ebounds_ref[0] = 0.01;
	     ebounds_ref[1] = 0.03;
	     ebounds_ref[2] = 0.07;
	     ebounds_ref[3] = 0.1;
	     ebounds_ref[4] = 0.3;
	     ebounds_ref[5] = 0.7;
	     ebounds_ref[6] = 1.0;
	     ebounds_ref[7] = 3.0;
	     ebounds_ref[8] = 7.0;
	     ebounds_ref[9] = 10.0;
	     ebounds_ref[10] = 30.0;
	     ebounds_ref[11] = 70.0;
	     ebounds_ref[12] = 100.0;


	     // Compare the grids.
	     if ( match( ebounds, ebounds_ref ) )
		 pass() << "Energy group boundary grid matches.";
	     else
		 fail() << "Energy group boundary grid did not match.";
	     
// 	     cout << "Energy Boundary Grid (keV):" << endl << endl;
// 	     for ( int i=0; i<ebounds.size(); ++i )
// 		 cout << "   ebounds[" << i << "] = " << ebounds[i] << endl;
	 }
     else
	 {
	     fail() << "The number of energy boundary points found in the data\n\t"
		    << "grid does not match the number returned by the\n\t"
		    << "getNumGroupBoundaries() accessor.";
	     fail() << "Did not test the results returned by\n\t"  
		    << "getGroupBoundaries().";
	 }
    
} // end of testEnergyBoundaryAccessor( SP<Opacity> spOpacity )

} // end namespace rtt_cdi_gandolf_test

//---------------------------------------------------------------------------//
//                            end of tCDIGandolf.cc
//---------------------------------------------------------------------------//
