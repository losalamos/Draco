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
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

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
    
    // Start the test.

//     cout << endl << "Testing the cdi_gandolf package." << endl;
//     cout << endl << "Create SP to a cdi_gandolf object." << endl << endl;

    SP<rtt_cdi_gandolf::GandolfFile> spGF;
    spGF = new rtt_cdi_gandolf::GandolfFile( op_data_file );

    if ( spGF )
	pass() << "SP to new GandolfFile object created for Al_BeCu.ipcress data.";
    else 
	{
	    fail() << "Failed to create SP to new GandolfFile object for Al_BeCu.ipcress data.";
	    return "Failed to finish testing cdi_gandolf.";
	}

    if ( spGF->getDataFilename() == op_data_file )
	pass() << "GandolfFile object is now linked to the Al_BeCu.ipcress data file.";
    else
	fail() << "GandolfFile object failed to link itself to theAl_BeCu.ipcress  data file.";

    if ( spGF->getNumMaterials() == 2 )
	pass() << "The correct number of materials was found in the Al_BeCu.ipcress data file.";
    else
	fail() << "spGF did not find the correct number of materials in the Al_BeCu.ipcress data file.";

    SP<rtt_cdi::Opacity> spOpacityABC;
    spOpacityABC = new rtt_cdi_gandolf::GandolfOpacity( spGF, matid );

    if ( spOpacityABC )
	pass() << "SP to new Opacity object created for Al_BeCu.ipcress data.";
    else
	fail() << "Failed to create SP to new Opacity object for Al_BeCu.ipcress data.";

    //------------------
    // Gray Opacity Test
    //------------------

    // Currently only Rossland opacities are available.  I need to
    // implement a mechanism so that the user can choose between
    // rosseland and plank.

    // Obtain a Rosseland Gray Opacity value for T=0.1 keV and density 
    // = 27.0 g/cm^3.
    double grayRosselandOpacity = spOpacityABC->getGrayRosseland( 0.1, 27.0 );

    // Print some info to standard out.
    // cout << endl << "grayRosselandOpacity for material Aluminum (matID=" 
    // << matid << ") is:" << endl
    // << "   opacity ( T=0.1 keV, rho=27.0 g/cm^3 ) = " <<
    // grayRosselandOpacity << " cm^2/gm." << endl << endl;


    // Make sure that the interpolated value matches previous
    // interpolations. 
    double tabulatedGrayOpacity = 4271.7041147070677; // cm^2/g

    if ( match( grayRosselandOpacity, tabulatedGrayOpacity ) )
	pass() << "grayRosselandOpacity computation was good for Al_BeCu.ipcress data.";
    else
	fail() << "grayRosselandOpacity value is out of spec for Al_BeCu.ipcress data.";
    
  //----------------
  // MG Opacity test
  //----------------
  
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
     
     // Interpolate the multigroup opacities for T = 0.01 keV and
     // density = 2.0 g/cm^3. 
     vector<double> mgRosselandOpacity
 	= spOpacityABC->getMGRosseland( 0.01, 2.0 );

//     // print the results
//     cout << endl << "Multigroup Rosseland Opacities for Aluminum (matID=" 
// 	 << matid << ") test:" << endl
// 	 << "   MGOpacities at T=0.01 keV and rho=2.0 g/cm^3 ( in cm^2/gm. )."
// 	 << endl << endl;
//     for ( int i=0; i<mgRosselandOpacity.size(); ++i )
//     cout << "   Opacity(group=" << std::setw(2) << i << ") = " 
// 	    << std::scientific << std::setprecision(5) 
// 	    << mgRosselandOpacity[i] << endl;
//     cout << endl;
    

//     // Compare the interpolated value with previous interpolations:
     if ( match( mgRosselandOpacity, tabulatedMGOpacity ) )
 	pass() << "Multigroup Opacity computation was good for Al_BeCu.ipcress data.";
     else
 	fail() << "Multigroup Opacity ocmputation failed for Al_BeCu.ipcress data.";

// -----------------------------------------------
// Test the data file "analyticOpacities.ipcress"
// -----------------------------------------------

     // The Opacities in this file are computed from the following
     // analytic formula:
     //     opacity = rho * T^4,
     // rho is the density and T is the temperature.

     // Gandolf data filename (IPCRESS format required)
     op_data_file = "analyticOpacities.ipcress";
     
     // Start the test.

     // SP<rtt_cdi_gandolf::GandolfFile> spGF;
     spGF = new rtt_cdi_gandolf::GandolfFile( op_data_file );
     
     if ( spGF )
	 pass() << "SP to new GandolfFile object created.";
     else
	 fail() << "Failed to create SP to new GandolfFile object.";
     
     if ( spGF->getDataFilename() == op_data_file )
	 pass() << "GandolfFile object is now linked to the data file.";
     else
	 fail() << "GandolfFile object failed to link itself to the data file.";

     if ( spGF->getNumMaterials() == 1 )
	 pass() << "The correct number of materials was found in the data file.";
     else
	 fail() << "spGF did not find the correct number of materials in the data file.";

     // SP<rtt_cdi::Opacity> spOpacityABC;
     spOpacityABC = new rtt_cdi_gandolf::GandolfOpacity( spGF, matid );
     
     if ( spOpacityABC )
	 pass() << "SP to new Opacity object created.";
     else
	 fail() << "Failed to create SP to new Opacity object.";

     //------------------
     // Gray Opacity Test
     //------------------
     
     // Currently only Rossland total opacities are available.  The
     // analytic formulat shown above does not include the addition of 
     // scattering.  I need to implement a mechanism so that the user
     // can choose absoptive-only Rosseland or Plank opacities.
     
     // Obtain a Rosseland Gray Opacity value for T=0.2 keV and density 
     // = 0.3 g/cm^3.
     double T = 1.0;   // keV
     double rho = 1.0; // g/cm^3
     grayRosselandOpacity = spOpacityABC->getGrayRosseland( T, rho );

    // Print some info to standard out.
//      cout << endl << "grayRosselandOpacity for material 1 (matID=" 
// 	  << matid << ") is:" << endl
// 	  << "   opacity ( T=" << T << " keV, rho=" << rho << " g/cm^3 ) = " <<
// 	 grayRosselandOpacity << " cm^2/gm." << endl << endl;

    // Make sure that the interpolated value matches previous
    // interpolations. 
     double scatter = 0.4;
     tabulatedGrayOpacity = rho*pow(T,4) + scatter; // cm^2/g

    if ( match( grayRosselandOpacity, tabulatedGrayOpacity ) )
	pass() << "grayRosselandOpacity computation was good for analyticOpacity data.";
    else
	fail() << "grayRosselandOpacity value is out of spec for analyticOpacity data.";
    
  //----------------
  // MG Opacity test
  //----------------
  
    // The solution to compare against:
    numGroups = 12;
    tabulatedMGOpacity.resize(numGroups);

    T   = 0.3; // keV
    rho = 0.7; // g/cm^3
    for ( int i=0; i<numGroups; ++i )
	tabulatedMGOpacity[i] = rho*pow(T,4); // cm^2/gm
    
    // Interpolate the multigroup opacities for T = 0.3 keV and
    // density = 0.7 g/cm^3. 
    mgRosselandOpacity
 	= spOpacityABC->getMGRosseland( T, rho );

//     // print the results
//     cout << endl << "Multigroup Rosseland Opacities for Aluminum (matID=" 
// 	 << matid << ") test:" << endl
// 	 << "   MGOpacities at T=0.01 keV and rho=2.0 g/cm^3 ( in cm^2/gm. )."
// 	 << endl << endl;
//     for ( int i=0; i<mgRosselandOpacity.size(); ++i )
//     cout << "   Opacity(group=" << std::setw(2) << i << ") = " 
// 	    << std::scientific << std::setprecision(5) 
// 	    << mgRosselandOpacity[i] << endl;
//     cout << endl;
    

//     // Compare the interpolated value with previous interpolations:
     if ( match( mgRosselandOpacity, tabulatedMGOpacity ) )
 	pass() << "Multigroup Opacity computation was good for analyticOpacity data.";
     else
 	fail() << "Multigroup Opacity ocmputation failed for analyticOpacity data.";


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

} // end namespace rtt_cdi_gandolf_test

//---------------------------------------------------------------------------//
//                            end of tCDIGandolf.cc
//---------------------------------------------------------------------------//
