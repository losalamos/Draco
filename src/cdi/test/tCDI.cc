//----------------------------------*-C++-*----------------------------------//
// tCDI.cc
// Kelly Thompson
// Thu Jun 22 13:07:00 2000
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// revision history:
// -----------------
// 1.1) Original

#include "tCDI.hh"
#include "../CDI.hh"
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
	using rtt_CDI_test::tCDI;
	return rtt_dsxx::SP<TestApp> ( new tCDI( argc, argv, os_in ));
    }
} // end namespace rtt_UnitTestFrame

// tCDI Stuff
//--------------------------------------------------
namespace rtt_CDI_test {

using std::string;
using std::cout;
using std::endl;
using std::vector;
using rtt_dsxx::SP;
using rtt_cdi::CDI;

tCDI::tCDI( int argc, char *argv[], std::ostream& os_in )
    : rtt_UnitTestFrame::TestApp( argc, argv, os_in )
{
    os() << "Created tCDI" << endl;
}

string tCDI::version() const
{
    return rtt_cdi::release();
}

string tCDI::runTest()
{
    
    // Get opacities from Gandolf
    // OpType must be one of { Gandolf, EOSPAC, Analytic }.
    rtt_cdi::OpType op_type = rtt_cdi::Gandolf;
    
    // Gandolf data filename (IPCRESS format required)
    string op_data_file = "Al_BeCu.ipcress";
    
    std::ifstream infile( op_data_file.c_str() );
    if ( ! infile )
	fail() << "Could not open file for reading.";
    else { 
	pass() << "File found for reading.";
	
	
	
	// Start the test.
	cout << endl
	     << "Testing the CDI package."
	     << endl;
	
	cout << endl
	     << "Create SP to a CDI object." << endl
	     << endl;
	
	SP<CDI> spCDI_Al;
	spCDI_Al = new CDI( op_type, op_data_file );
	
	cout << endl
	     << "The Gandalf input file is named: \"" 
	     << spCDI_Al->getOpacityDataFilename() << "\""
	     << endl << endl;
	
	vector<int> matids = spCDI_Al->getMatIDs();
	
	cout << endl << "Material IDs found:" << endl;
	for ( int i=0; i < matids.size(); ++i )
	    cout << "   " << matids[i] << endl;
	
	
    
    
	
    //    double grayRosselandOpacity = spCDI_Al->getGrayOpacity( 1.0, 10.0 );

	pass() << "Done testing CDI.";

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

} // end namespace rtt_CDI_test

//---------------------------------------------------------------------------//
//                            end of tCDI.cc
//---------------------------------------------------------------------------//
