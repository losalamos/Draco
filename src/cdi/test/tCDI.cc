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
#include "DummyOpacity.hh"
#include "../CDI.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

#include <iostream>
#include <vector>
//#include <fstream>

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
    
    // Start the test.
    cout << endl
	 << "Testing the CDI package."
	 << endl;
	
    cout << endl
	 << "Create SP to a Opacity object." << endl
	 << endl;
    
    SP<rtt_cdi::Opacity> spOpacity;
    spOpacity = new rtt_dummy_opacity::DummyOpacity();

    cout << endl
	 << "Create SP to a CDI object." << endl
	 << endl;
    
    SP<CDI> spCDI_mat1;
    spCDI_mat1 = new CDI( spOpacity );
	
    // do some dummy calls here to make sure things are working
    //
    // vector<double> opacities 
    //    = spCDI_mat1->getGrayOpacity( double temp, 
    //                                  double density);

	pass() << "Done testing CDI.";

	cout << endl << endl;


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
