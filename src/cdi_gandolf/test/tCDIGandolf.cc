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
    
    // Get opacities from Gandolf
    // OpType must be one of { Gandolf, EOSPAC, Analytic }.
    // rtt_cdi::OpType op_type = rtt_cdi::Gandolf;
    
    // Gandolf data filename (IPCRESS format required)
    // string op_data_file = "../../../../../gandolf/ipcress/Al_BeCu.ipcress";
    string op_data_file = "Al_BeCu.ipcress";
    
    std::ifstream infile( op_data_file.c_str() );
    if ( ! infile )
	fail() << "Could not open file for reading.";
    else { 
	pass() << "File found for reading.";
	
	
	
	// Start the test.
	cout << endl
	     << "Testing the cdi_gandolf package."
	     << endl;
	
	cout << endl
	     << "Create SP to a cdi_gandolf object." << endl
	     << endl;
	
	//SP<CDI> spCDI_ABC;
	//spCDI_ABC = new CDI( op_type, op_data_file );
	
	SP<rtt_cdi::Opacity> spOpacityABC;
	spOpacityABC = new rtt_cdi_gandolf::GandolfOpacity( op_data_file );

	cout << endl
	     << "The Gandalf input file is named: \"" 
	     << spOpacityABC->getDataFilename() << "\""
	     << endl << endl;
	
	vector<int> matids = spOpacityABC->getMatIDs();
	
	cout << endl << "Material IDs found:" << endl;
	for ( int i=0; i < matids.size(); ++i )
	    cout << "   " << matids[i] << endl;
	
	
    
    
	
    //    double grayRosselandOpacity = spCDI_ABC->getGrayOpacity( 1.0, 10.0 );

	pass() << "Done testing cdi_gandolf.";

	cout << endl << endl;
    }


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
