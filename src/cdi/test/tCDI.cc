//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/tCDI.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 13:07:00 2000
 * \brief  Implementation file for the CDI unit test.
 */
//---------------------------------------------------------------------------//
// $Id$
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
	
    // Create an opacity object
    SP<rtt_cdi::Opacity> spOpacity;
    if ( spOpacity = new rtt_dummy_opacity::DummyOpacity() )
	pass() << "SP to opacity object created successfully.";
    else
	fail() << "Failed to create SP to opacity object.";

    // Create a CDI object linked to the opacity object    
    SP<CDI> spCDI_mat1;
    if ( spCDI_mat1 = new CDI( spOpacity ) )
	pass() << "SP to CDI object created successfully.";
    else
	fail() << "Failed to create SP to CDI object.";
	

    // test the getDataFilename() function.
    string fname = spOpacity->getDataFilename();
    cout << "The data file is named: " << fname << endl;
    fname.append("blah");
    cout << "The data file is named: " << fname << endl;



    //----------------------------------------
    // Start the tests
    //----------------------------------------

    // do some dummy calls here to make sure things are working
    // These values are actually ignored by DummyOpacity.

    double temp = 1.0;        // keV
    double density = 27.0;    // g/cm^3

    // --> Try to collect gray opacity data.

    double grayOpacityReference = temp + density/10000;  // cm^2/g
    double grayOpacity 
        = spCDI_mat1->getGrayRosselandOpacity( temp, density);

    if ( match( grayOpacity, grayOpacityReference ) )
	pass() << "Access to gray opacity data succeeded.";
    else
	fail() << "Access to gray opacity data failed.";

    // --> Try to collect multigroup opacity data.

    // In DummyOpacity ngroups is hard coded to 3.
    vector<double> MGOpacitiesReference(3); 
    for ( int i=0; i<3; ++i )
	MGOpacitiesReference[i] = (i+1)*1000.0 + temp + density/10000;

    vector<double> MGOpacities 
	= spCDI_mat1->getMGRosselandOpacity( temp, density);

    if ( match( MGOpacities, MGOpacitiesReference ) )
	pass() << "Access to multigroup opacity data succeeded.";
    else
	fail() << "Access to multigroup data failed.";

    //----------------------------------------
    // End of tests
    //----------------------------------------
    
    pass() << "Done testing CDI.";
    cout << endl << endl;

    //----------------------------------------
    // Print the test result.
    //----------------------------------------

    if (passed()) {
	pass() << "All tests passed.";
	return "All tests passed.";
    }
    return "Some tests failed.";
}

//---------------------------------------------
// Compare Reference value to computed values
//---------------------------------------------
bool tCDI::match( const vector<double> computedValue, 
			 const vector<double> referenceValue )
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
} // end of tCDI::match( vector<double>, vector<double> )

bool tCDI::match( const double computedValue,
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

} // end of tCDI::match( double, double )

} // end namespace rtt_CDI_test

//---------------------------------------------------------------------------//
//                            end of tCDI.cc
//---------------------------------------------------------------------------//
