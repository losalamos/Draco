//----------------------------------*-C++-*----------------------------------//
/*!
 * \file tCDIGandolf.cc
 * \author Kelly Thompson
 * \date Thu Jun 22 13:07:00 2000
 * \brief 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "tCDIGandolf.hh"
#include "../GandolfOpacity.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>

// Unit Test Frame Stuff
//----------------------------------------
namespace rtt_UnitTestFrame {
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
					   std::ostream& os_in ) {
	using rtt_cdi_gandolf_test::tCDIGandolf;
	return rtt_dsxx::SP<TestApp> ( 
	    new tCDIGandolf( argc, argv, os_in ));
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
//using rtt_cdi::Opacity;
//using rtt_cdi::CDI;

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

    const double TOL = 1.0e-10;
    
    // Get opacities from Gandolf
    // OpType must be one of { Gandolf, EOSPAC, Analytic }.
    // rtt_cdi::OpType op_type = rtt_cdi::Gandolf;
    
    // Gandolf data filename (IPCRESS format required)
    // string op_data_file = "../../../../../gandolf/ipcress/Al_BeCu.ipcress";
    string op_data_file = "Al_BeCu.ipcress";
    int matid=10001;
    
//     std::ifstream infile( op_data_file.c_str() );
//     if ( ! infile )
// 	fail() << "Could not open file for reading.";
//     else { 
// 	pass() << "File found for reading.";
	
	// Start the test.
	cout << endl
	     << "Testing the cdi_gandolf package."
	     << endl;
	
	cout << endl
	     << "Create SP to a cdi_gandolf object." << endl
	     << endl;
	
	SP<rtt_cdi::Opacity> spOpacityABC;
	spOpacityABC = new 
	    rtt_cdi_gandolf::GandolfOpacity( op_data_file, matid );

	cout << endl
	     << "The Gandalf input file is named: \"" 
	     << spOpacityABC->getDataFilename() << "\""
	     << endl;

	cout << "Material " << matid << " was found in the data file."
	     << endl << endl;


	// how to choose between rosseland and plank?
	double grayRosselandOpacity = spOpacityABC->getGray( 0.1, 27.0 );

	cout << endl << "grayRosselandOpacity for material Aluminum (matID=" 
	     << matid << ") is:" << endl
	     << "   opacity ( T=0.1 keV, rho=27.0 g/cm^3 ) = " <<
	    grayRosselandOpacity << " cm^2/gm." << endl << endl;

	double tabulatedGrayOpacity = 4271.7041147070677; // cm^2/g
	double absdiff = fabs( ( grayRosselandOpacity
				 - tabulatedGrayOpacity )
			       / tabulatedGrayOpacity );

	if ( absdiff > TOL )
	    fail() << "grayRosselandOpacity value is out of spec.";
	else
	    pass() << "grayRosselandOpacity computation was good.";
	
	// MG Opacity test.

	cout << endl << "Multigroup Rosseland Opacities for Aluminum (matID=" 
	     << matid << ") test:" << endl
	     << "   MGOpacities at T=0.01 keV and rho=2.0 g/cm^3 ( in cm^2/gm. )."
	     << endl << endl;

	// The solution to compare against:
	const int numGroups = 33;
	const double tabulatedMGOpacity[] =
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
	
	vector<double> mgRosselandOpacity
	    = spOpacityABC->getMG( 0.01, 2.0 );

// 	cout << endl;
// 	cout << "Group, mgRosselandOpacity, tabulatedMGOpacity, absdiff" << endl;

	bool exactMatch = true;
	for ( int i=0; i<numGroups; ++i ) 
	    {
		absdiff = fabs ( ( mgRosselandOpacity[i] 
				   - tabulatedMGOpacity[i] )
				 / tabulatedMGOpacity[i] );
		
 		cout << "   MGOpacity(group="<< std::setw(2) << i+1 << ") = " 
 		     << std::scientific << std::setprecision(5) 
 		     << mgRosselandOpacity[i] << " cm^2/gm." << endl;

		if ( absdiff > TOL ) 
		    {
			exactMatch = false;
			break;
		    }
	    }

	if ( exactMatch )
	    {
		pass() << "Multigroup Opacity computation was good.";
	    }
	else
	    fail() << "Multigroup Opacity ocmputation failed.";

	cout << endl << endl;

	// } // data file found


    // Print the test result.
    // ----------------------------------------

    if (passed()) {
	pass() << "All tests passed.";
	return "All tests passed.";
    }
    return "Some tests failed.";
}

} // end namespace rtt_cdi_gandolf_test

//---------------------------------------------------------------------------//
//                            end of tCDIGandolf.cc
//---------------------------------------------------------------------------//
