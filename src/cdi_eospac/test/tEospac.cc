//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_eospac/test/tEospac.cc
 * \author Kelly Thompson
 * \date   Mon Apr 2 14:20:14 2001
 * \brief  Implementation file for tEospac
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "tEospac.hh"
//#include "tEospac.t.hh"

#include "../Eospac.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

//#include <vector>

// for debuging only.
//#include <iomanip>
#include <iostream>

// Unit Test Frame Stuff
//----------------------------------------
namespace rtt_UnitTestFrame {
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
					   std::ostream& os_in ) {
	return rtt_dsxx::SP<TestApp> ( 
	    new rtt_cdi_eospac_test::tEospac( argc, argv, os_in ));
    }
} // end namespace rtt_UnitTestFrame

// tEospac Stuff
//--------------------------------------------------
namespace rtt_cdi_eospac_test {
    
    tEospac::tEospac( int argc, char *argv[], std::ostream& os_in )
	: rtt_UnitTestFrame::TestApp( argc, argv, os_in )
	{
	    os() << "Created tEospac" << std::endl;
	}
    
    std::string tEospac::version() const
	{
	    return rtt_cdi_eospac::release();
	}
    
    //===========================================================================
    /*!
     * \brief
     *
     * On the XDIV LAN the EOSPAC library is located at:
     *
     * /usr/local/codes/data/eos/eospac_5-30beta/lib/sgi/64bit/libeospac.a
     */
    //===========================================================================
    std::string tEospac::runTest()
	{
	    
	    std::cout << std::endl
		 << "Test of C++ code calling EOSPAC routines" 
		 << std::endl << std::endl;

	    // Set the material identifier
	    // This one is for Aluminum (03717) 
	    // Category-1 data (0) + Mat# 371 (Al) + Version # 7

	    // See http://int.lanl.gov/projects/sdm/win/materials/ for 
	    // material ID information.

	    const int matID = 3717;

	    // Create an Eospac object
	    
	    rtt_dsxx::SP< rtt_cdi_eospac::Eospac > spEospac;
	    
	    if ( spEospac = new rtt_cdi_eospac::Eospac( matID ) )
		pass() << "SP to new Eospac object created.";
	    else
		{
		    fail() << "Unable to create SP to new Eospac object.";
		    return "Unable to create SP to new Eospac object.";
		}
	    
	    // Get an Electron internal energy value;

	    fail() << "change this to Cv for e-";

	    double density = 1.0; // (Mg/m^3)
	    double temperature = 5800; // K

	    double specificElectronInternalEnergy =
		spEospac->getSpecificElectronInternalEnergy(
		    density, temperature );

	    std::cout << "specificElectronInternalEnergy = " 
		      << specificElectronInternalEnergy << std::endl;

	    std::cout <<  std::endl <<  std::endl;

	    // ---------------------- //
	    // Print the test result. //
	    // ---------------------- //
	    
	    if ( passed() ) {
		pass() << "All tests passed.";
		return "All tests passed.";
	    }
	    
	    return "Some tests failed.";
	    
	} // end of runTest()
    
//     // ------------------------------------------ //
//     // Compare Reference value to computed values //
//     // ------------------------------------------ //
    
//     bool tEospac::match( const double computedValue,
// 			 const double referenceValue ) const
// 	{
// 	    // Start by assuming that the two quantities match exactly.
// 	    bool em = true;
	    
// 	    // Compare items up to 10 digits of accuracy.
	    
// 	    const double TOL = 1.0e-10;
	    
// 	    // Calculate the absolute value of the relative difference between 
// 	    // the computed and reference values.
	    
// 	    std::cout.precision(10);
// 	    std::cout << "\t" << computedValue
// 		      << "\t" << referenceValue << std::endl;
	    
// 	    double reldiff = fabs( ( computedValue - referenceValue )
// 				   / referenceValue );
	    
// 	    // If the comparison fails then change the value of "em" return
// 	    // the result;
// 	    if ( reldiff > TOL )
// 		em = false;
	    
// 	    return em;    
	    
// 	} // end of tEospac::match( double, double )
    
//     // ------------- //
//     // Match vectors //
//     // ------------- //
    
//     bool tEospac::match( 
// 	const std::vector< double >& computedValue, 
// 	const std::vector< double >& referenceValue ) const
// 	{
// 	    // Start by assuming that the two quantities match exactly.
// 	    bool em = true;
	    
// 	    // Compare items up to 10 digits of accuracy.
// 	    const double TOL = 1.0e-10;
	    
// 	    // Test each item in the list
// 	    double reldiff = 0.0;
// 	    for ( int i=0; i<computedValue.size(); ++i )
// 		{
// 		    std::cout.precision(10);
// 		    std::cout << "\t" << computedValue[i]
// 			      << "\t" << referenceValue[i] << std::endl;
		    
// 		    reldiff = fabs( ( computedValue[i] - referenceValue[i] )
// 				    / referenceValue[i] );
// 		    // If the comparison fails then change the value of "em"
// 		    // and exit the loop.
		    
// 		    // DEBUG: must #include <iomanip>
// 		    //
// 		    // 	    std::cout << std::setprecision(14) << "   "
// 		    // 		      << computedValue[i] << "   "
// 		    // 		      << referenceValue[i] << "   "
// 		    // 		      << reldiff << std::endl;
		    
// 		    if ( reldiff > TOL )
// 			{
// 			    em = false;
// 			    break;
// 			}
// 		}
// 	    return em;
// 	} 
    
//     // ------------------------ //
//     // Match vectors of vectors //
//     // ------------------------ //
    
//     bool tEospac::match( 
// 	const std::vector< std::vector<double> >& computedValue, 
// 	const std::vector< std::vector<double> >& referenceValue ) const 
// 	{
// 	    // Start by assuming that the two quantities match exactly.
// 	    bool em = true;
	    
// 	    // Compare items up to 10 digits of accuracy.
// 	    const double TOL = 1.0e-10;
	    
// 	    // Test each item in the list
// 	    double reldiff = 0.0;
// 	    for ( int i=0; i<computedValue.size(); ++i )
// 		{
// 		    for ( int j=0; j<computedValue[i].size(); ++j )
// 			{

// 			    std::cout.precision(10);
// 			    std::cout << "\t" << computedValue[i][j]
// 				      << "\t" << referenceValue[i][j] << std::endl;
			    
// 			    reldiff = fabs( ( computedValue[i][j] - referenceValue[i][j] )
// 					    / referenceValue[i][j] );
// 			    // If the comparison fails then change the value of "em"
// 			    // and exit the loop.
// 			    if ( reldiff > TOL ) {
// 				em = false; break; }
// 			}
// 		}
// 	    return em;
// 	} 
    
} // end namespace rtt_cdi_eospac_test

//---------------------------------------------------------------------------//
// end of tEospac.cc
//---------------------------------------------------------------------------//
