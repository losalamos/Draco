//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tDummyEoS.cc
 * \author Kelly Thompson
 * \date   Mon Apr 16 13:28:29 2001
 * \brief  Implementation file for tDummyEoS
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <vector>  // required by match()
#include <cmath>   // define fabs()

#include "tDummyEoS.hh"
#include "../Release.hh"

#include "DummyEoS.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

// --------------------- //
// Unit Test Frame Stuff //
// --------------------- //

namespace rtt_UnitTestFrame {
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
					   std::ostream& os_in ) {
	return rtt_dsxx::SP<TestApp> ( 
	    new rtt_dummy_eos_test::tDummyEoS( argc, argv, os_in ));
    }
} // end namespace rtt_UnitTestFrame


// ------------------- //
// tDummyEoS Stuff //
// ------------------- //

namespace rtt_dummy_eos_test 
{
    
    tDummyEoS::tDummyEoS( int argc, char *argv[], std::ostream& os_in )
	: rtt_UnitTestFrame::TestApp( argc, argv, os_in )
	{
	    os() << "Created tDummyEoS" << std::endl;
	}
    
    std::string tDummyEoS::version() const
	{
	    return rtt_cdi::release();
	}
    
    
    //=============================== //
    // The DummyEoS tests start here. //
    //=============================== //


    std::string tDummyEoS::runTest()
	{
	    
	    // --------------------- //
	    // Create an EoS object. //
	    // --------------------- //
	    
	    // The smart pointer points to a generic EoS object.
	    rtt_dsxx::SP< rtt_cdi::EoS > spEoS;
	    
	    // The actual instatniate is specific (dummyEoS).
	    if ( spEoS = new rtt_dummyEoS::DummyEoS() )
		// If we get here then the object was successfully instantiated.
		pass() << "Smart Pointer to new EoS object created.";
	    else
		{
		    fail() << "Unable to create a Smart Pointer to new EoS object.";
		    return "Unable to create a Smart Pointer to new EoS object.";
		}
	    
	    // --------- //
	    // EoS Tests //
	    // --------- //
	    
	    double temperature = 5800.0; // Kelvin
	    double density     = 27.0;   // g/cm^3
	    double tabulatedSpecificElectronInternalEnergy
		= temperature + 1000.0*density; // kJ/g
	    
	    double seie = spEoS->getSpecificElectronInternalEnergy( 
		temperature, density );
	    
	    if ( match ( seie, tabulatedSpecificElectronInternalEnergy ) ) 
// 		pass() << spEoS->getDataDescriptor()
		pass() << "The getSpecificElectronInternalEnergy( dbl, dbl)"
		       << " request returned the expected value.";
	    else
// 		fail() << spEoS->getDataDescriptor()
		fail() << "The getSpecificElectronInternalEnergy( dbl, dbl)"
		       << " request returned a value that is out of spec.";
	    
	    // try using a vectors of temps. and densities
	    // vtemperature.size() == vdensity.size()
	    
	    std::vector< double > vtemperature(3);
	    vtemperature[0] = 5000.0; // Kelvin
	    vtemperature[1] = 7000.0; // Kelvin
	    vtemperature[2] = 3000.0; // Kelvin

	    std::vector< double > vdensity(3);
	    vdensity[0] = 0.35; // g/cm^3
	    vdensity[1] = 1.0;  // g/cm^3
	    vdensity[2] = 9.8;  // g/mcm^3

	    // Retrieve electron based heat capacities.
	    std::vector< double > vRefCve( vtemperature.size() );
	    for ( int i=0; i<vtemperature.size(); ++i )
		vRefCve[i] = vtemperature[i] + vdensity[i]/1000.0;
	    
	    std::vector< double > vCve = spEoS->getElectronHeatCapacity(
		vtemperature, vdensity );
	    
	    if ( match ( vCve, vRefCve ) ) 
// 		pass() << spEoS->getDataDescriptor()
		pass() << "The getElectronHeatCapacity( vect, vect ) request"
		       << " returned the expected values.";
	    else
// 		fail() << spEoS->getDataDescriptor()
		fail() << "The getElectronHeatCapacity( vect, vect ) request"
		       << " returned values that are out of spec.";


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
    
    bool tDummyEoS::match( double computedValue,
			   double referenceValue ) const
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
	    
	} // end of tDummyEoS::match( double, double )
    
    bool tDummyEoS::match( const std::vector< double >& computedValue, 
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
	} // end of tDummyEoS::match( vector<double>, vector<double> )
    
} // end namespace rtt_dummy_eos_test

//---------------------------------------------------------------------------//
// end of tDummyEoS.cc
//---------------------------------------------------------------------------//
