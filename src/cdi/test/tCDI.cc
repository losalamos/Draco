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
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

#include "DummyGrayOpacity.hh"
#include "DummyMultigroupOpacity.hh"
#include "../CDI.hh"


// --------------------- //
// Unit Test Frame Stuff //
// --------------------- //

namespace rtt_UnitTestFrame {
    rtt_dsxx::SP< TestApp > 
	TestApp::create( int &argc, char *argv[],
			 std::ostream& os_in ) 
	{
	    return rtt_dsxx::SP< TestApp > 
		( new rtt_CDI_test::tCDI( argc, argv, os_in ) );
	}
} // end namespace rtt_UnitTestFrame

// ---------- //
// tCDI Stuff //
// ---------- //

namespace rtt_CDI_test {

tCDI::tCDI( int argc, char *argv[], std::ostream& os_in )
    : rtt_UnitTestFrame::TestApp( argc, argv, os_in )
{
    os() << "Created tCDI" << std::endl;
}

std::string tCDI::version() const
{
    return rtt_cdi::release();
}



		 // ============================== //
		 // The CDI unit tests start here. //
		 // ============================== //


std::string tCDI::runTest()
{
    // Start the test.

    std::cout << std::endl
	      << "Testing the CDI package."
	      << std::endl;
	

    // ----------------------------------- //
    // Create some Opacity and CDI objects //
    // ----------------------------------- //

    // Create a GrayOpacity object

    rtt_dsxx::SP< rtt_cdi::GrayOpacity > spDGrO;
    
    if ( spDGrO = new rtt_dummyGrayOpacity::DummyGrayOpacity() )
	pass() << "SP to new GrayOpacity object created.";
    else
	{
	    fail() << "Unable to create a SP to new GrayOpacity object.";
	    return "Unable to create a SP to new GrayOpacity object.";
	}


    // Create a MultigroupOpacity object

    rtt_dsxx::SP< rtt_cdi::MultigroupOpacity > spDMgO;
    
    if ( spDMgO = new rtt_dummyMultigroupOpacity::DummyMultigroupOpacity() )
	pass() << "SP to new MultigroupOpacity object created.";
    else
	{
	    fail() << "Unable to create a SP to new MultigroupOpacity object.";
	    return "Unable to create a SP to new MultigroupOpacity object.";
	}

    // Create a CDI object linked to a DummyGrayOpacity object.

     rtt_dsxx::SP< rtt_cdi::CDI > spCdiDumGr;
     if ( spCdiDumGr = new rtt_cdi::CDI( spDGrO ) )
 	pass() << "SP to CDI object created successfully (GrayOpacity).";
     else
	 fail() << "Failed to create SP to CDI object (GrayOpacity).";
     
     // Create a CDI object linked to a DummyMultigroupOpacity object.

     rtt_dsxx::SP< rtt_cdi::CDI > spCdiDumMg;
     if ( spCdiDumMg = new rtt_cdi::CDI( spDMgO ) )
 	pass() << "SP to CDI object created successfully (MultigroupOpacity).";
     else
	 fail() << "Failed to create SP to CDI object (MultigroupOpacity).";


     // --------------- //
     // Start the tests // 
     // --------------- //

     // gray test case:  Find the value of opacity at T=0.35 keV and
     //      rho = 27.2 g/cm^3.  For DummyGrayOpacity the value
     //      should be .35272 cm^2/g.
     
     double temp = 0.35; // keV
     double dens = 27.2; // g/cm^3
     double refOpacity = temp + dens/1000.0; // cm^2/g

     double opacity = spCdiDumGr->gray()->getOpacity( temp, dens );

     if ( match( opacity, refOpacity ) )
	 pass() << "CDI->gray()->getOpacity(T,rho) computation was good.";
     else
	 fail() << "CDI->gray()->getOpacity(T,rho) value is out of spec.";


     // mg test case:  Find the mg opacities at T=0.35 keV and 
     //     rho = 27.2 g/cm^3.  For DummyMultigroupOpacity the values
     //     should be { }
     
     int ng = spCdiDumMg->mg()->getNumGroups();

     if ( ng == 3 ) // DummyMultigroupOpacity is hardwired to 3 groups.
	 pass() << "CDI->mg()->getNumGroups() access was good.";
     else
	 fail() << "CDI->mg()->getNumGroups() access failed.";

     // The energy groups in DummyMultigroupOpacity are hardwired to
     // be { 0.05, 0.5, 5.0, 50.0 } keV.
     std::vector< double > refEnergyBoundary(ng+1);
     refEnergyBoundary[0] = 0.05;
     refEnergyBoundary[1] = 0.5;
     refEnergyBoundary[2] = 5.0;
     refEnergyBoundary[3] = 50.0;

     std::vector< double > energyBoundary(ng+1);
     energyBoundary = spCdiDumMg->mg()->getGroupBoundaries();

     if ( match( energyBoundary, refEnergyBoundary ) )
	 pass() << "Accessed energy group boundaries for\n"
		<< "\tCDI->mg()->getGroupBoundaries() match expected values.";
     else
	 fail() << "Accessed energy group boundaries for\n"
		<< "\tCDI->mg()->getGroupBoundaries() do not match expected values.";

     std::vector< double > vRefOpacity( ng );
     for ( int ig=0; ig<ng; ++ig )
	 {
	     vRefOpacity[ig] = 2.0*(temp+dens/1000.0)
		 /(energyBoundary[ig]+energyBoundary[ig+1]);
	     std::cout << "vRefOpacity[" << ig << "] = "
		       << vRefOpacity[ig] << " cm^2/g." << std::endl;
	 }

     std::vector< double > vOpacity( ng );
     vOpacity = spCdiDumMg->mg()->getOpacity( temp, dens );

     if ( match( vOpacity, vRefOpacity ) )
	 pass() << "CDI->mg()->getOpacity(T,rho) computation was good.";
     else
	 fail() << "CDI->mg()->getOpacity(T,rho) computation is out of spec.";


    //----------------------------------------
    // End of tests
    //----------------------------------------
    
    pass() << "Done testing CDI.";
    std::cout << std::endl << std::endl;

    //----------------------------------------
    // Print the test result.
    //----------------------------------------

    if ( passed() ) {
	pass() << "All tests passed.";
	return "All tests passed.";
    }
    return "Some tests failed.";
}

bool tCDI::match( const double computedValue,
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

} // end of tCDI::match( double, double )

bool tCDI::match( const std::vector< double >& computedValue, 
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
} // end of tCDI::match( vector<double>, vector<double> )

} // end namespace rtt_CDI_test

//---------------------------------------------------------------------------//
//                            end of tCDI.cc
//---------------------------------------------------------------------------//
