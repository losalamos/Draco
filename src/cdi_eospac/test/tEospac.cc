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

// cdi_gandolf dependencies
#include "tEospac.hh"
#include "../Eospac.hh"
#include "../SesameTables.hh"
#include "../Release.hh"

// Draco dependencies
#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

// STL dependencies
#include <iostream>
#include <cmath>     // define fabs()

// DEBUG dependencies
// #include <iomanip>

// Unit Test Frame Stuff
//----------------------------------------
namespace rtt_UnitTestFrame 
{
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
					   std::ostream& os_in )
	{
	    return rtt_dsxx::SP<TestApp> ( 
		new rtt_cdi_eospac_test::tEospac( argc, argv, os_in ) );
	}
} // end namespace rtt_UnitTestFrame

// tEospac Stuff
//--------------------------------------------------
namespace rtt_cdi_eospac_test 
{
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
     * \brief Tests the Eospac constructor and access routines.
     *
     * On the XDIV LAN the EOSPAC library is located at:
     *
     * /usr/local/codes/data/eos/eospac_5-30beta/lib/sgi/64bit/libeospac.a
     *
     * We have a slightly modified copy (added one routine to help
     * C/F77 translation of character arrays) located at:
     *
     * /radtran/vendors/eospac/IRIX64/lib64/libeospac.a
     *
     * To use this package Draco must be compiled with the 
     * --with-eospac-lib=/radtran/vendors/eospac/IRIX64/lib64 tag.
     *
     */
    //===========================================================================
    std::string tEospac::runTest()
	{
	    
	    std::cout << std::endl
		 << "Test of C++ code calling EOSPAC routines" 
		 << std::endl << std::endl;


	    // ---------------------------- //
	    // Create a SesameTables object //
	    // ---------------------------- //

	    // The user must create a SesameTables object that links
	    // material ID numbers to EOSPAC data types (each
	    // SesameTables object only contains lookups for one
	    // material).   If the user needs heat capacity values
	    // for Al then he/she must create a SesameTables object
	    // for Aluminum and then assign an aluminum material ID
	    // (e.g. 3717) to the enelc EOSPAC data type.  See the
	    // tests below for more details. 

	    // Set the material identifier
	    // This one is for Aluminum (03717) 
	    // Category-1 data (0) + Mat# 371 (Al) + Version # 7

	    // See http://int.lanl.gov/projects/sdm/win/materials/ for 
	    // material ID information.

	    // This matID for Al has lookup tables for prtot, entot,
	    // tptot, tntot, pntot, eptot, prelc, enelc, tpelc, tnelc
	    // pnelc, epelc, prcld, and encld (see SesameTables.hh for 
	    // an explanantion of these keywords).  I need the table
	    // that contains enion lookups so that I can query for
	    // Cve() values.

	    const int Al3717 = 3717;

	    // This matId for Al has tables 101, 102, 201, 601, 602,
	    // 603 and 604.  I need table 601 for zfree.

	    const int Al23714 = 23714;

	    // Create a SesameTables object for Aluminum.

	    rtt_cdi_eospac::SesameTables AlSt;

	    // Assign matID Al3717 to enion lookups (used for Cvi) for 
	    // AlSt.  We can also assign these tables when the Eospac
	    // object is created (see example below). 

	    // Also assign matID Al23714 for temperature-based
	    // electron thermal conductivity (tconde).

	    AlSt.enion( Al3717 ).tconde( Al23714 );

	    // Verify that the assignments were made correctly.

	    // Cvi (returnType=8=ES4enion) should point to matID
	    // 3717.  The user should never need to access this
	    // function.  However Eospac.cc does and we need to test
	    // this funcitonality.

	    if ( AlSt.matID( rtt_cdi_eospac::ES4enion ) != 3717 )
		fail() << "AlSt.matID(ES4enion) points to the wrong matID.";

	    // The temperature-based electorn thermal conductivity
	    // (returnType=27=ES4tconde) should point to matID
	    // 23714.  The user should never need to access this
	    // function.  However Eospac.cc does and we need to test
	    // this funcitonality.

	    if ( AlSt.matID( rtt_cdi_eospac::ES4tconde ) != 23714 )
		fail() << "AlSt.matID(27) points to the wrong matID.";	    
	    

	    // ----------------------- //
	    // Create an Eospac object //
	    // ----------------------- //
	    

	    // An Eospac object allows the user to access EoS
	    // information about a material that has been constructed 
	    // in a SesameTable object.  The constructor for Eospac
	    // takes one argument: a SesameTables object.
	    
	    rtt_dsxx::SP< const rtt_cdi_eospac::Eospac > spEospac;

	    // Try to instantiate the new Eospac object.
	    // Simultaneously, we are assigned material IDs to more
	    // SesameTable values.
	    
	    if ( spEospac = new rtt_cdi_eospac::Eospac( 
		    AlSt.enelc( Al3717 ).zfree3( Al23714 ) ) )

		// Alternatively, we can avoid carrying around the
		// AlSt object.  We can, instead, create a temporary
		// version that is only used here in the constructor
		// of Eospac().		

//		spEospac = new rtt_cdi_eospac::Eospac( 
// 		rtt_cdi_eospac::SesameTables().enelc( Al3717 )
//                  .zfree3( Al23714 ).enion( Al3717 ).tconde( Al23714 ) ) )

		pass() << "SP to new Eospac object created.";
	    else
		{
		    fail() << "Unable to create SP to new Eospac object.";
		    return "Unable to create SP to new Eospac object.";
		}
	    
	    // --------------------------- //
	    // Test scalar access routines //
	    // --------------------------- //

	    const double K2keV = 1.0/1.1604412E+7; // keV/Kelvin

	    // All of these tests request an EoS value given a single
	    // temperature and a single density.

	    // Retrieve an Electron internal energy value;

	    double density     = 1.0;  // g/cm^3
	    double temperature = 5800; // K
	    temperature *= K2keV;      // convert temps to keV

	    double refValue = 1.052552479800656;  // kJ/g

	    double specificElectronInternalEnergy =
		spEospac->getSpecificElectronInternalEnergy(
		    temperature, density );

	    if ( match( specificElectronInternalEnergy, refValue ) )
		pass() << "getSpecificElectronInternalEnergy() test passed.";
	    else
		fail() << "getSpecificElectronInternalEnergy() test failed.";

	    // Retrieve an electron heat capacity (= dE/dT)	    

	    refValue = 3146.719924188898; // kJ/g/keV
	    
	    double heatCapacity =
		spEospac->getElectronHeatCapacity( temperature,
						   density );

	    if ( match(  heatCapacity, refValue ) )
		pass() << "getElectronHeatCapacity() test passed.";
	    else
		fail() << "getElectronHeatCapacity() test failed.";

	    // Retrive an Ion Internal Energy

	    refValue = 5.238217222081386; // kJ/g

	    double specificIonInternalEnergy = 
		spEospac->getSpecificIonInternalEnergy( 
		    temperature, density );

	    if ( match( specificIonInternalEnergy, refValue ) )
		pass() << "getSpecificIonInternalEnergy() test passed.";
	    else
		fail() << "getSpecificIonInternalEnergy() test failed.";

	    // Retrieve an ion based heat capacity

	    refValue = 6748.931926862662; // kJ/g/keV

	    heatCapacity =
		spEospac->getIonHeatCapacity( temperature, density );
	    
	    if ( match( heatCapacity, refValue ) )
		pass() << "getIonHeatCapacity() test passed.";
	    else
		fail() << "getIonHeatCapacity() test failed.";

	    // Retrieve the number of free electrons per ion

	    refValue = 12.89854626207534; // electrons per ion

	    double nfree =
		spEospac->getNumFreeElectronsPerIon( 
		    temperature, density );
	    
	    if ( match( nfree, refValue ) )
		pass() << "getNumFreeElectronsPerIon() test passed.";
	    else
		fail() << "getNumFreeElectronsPerIon() test failed.";

	    // Retrieve the electron based thermal conductivity

	    refValue = 1.389598060091371e+29; // 1/s/cm

	    double chie = 
		spEospac->getElectronThermalConductivity(
		    temperature, density );

	    if ( match( chie, refValue ) )
		pass() << "getElectronThermalConductivity() test passed.";
	    else
		fail() << "getElectronThermalConductivity() test failed.";

	    // --------------------------- //
	    // Test vector access routines //
	    // --------------------------- //

	    // Set up simple temp and density vectors.  vtemp(i) will
	    // always be associated with vdensities(i).  In this case
	    // both tuples have identical data so that the returned
	    // results will also be identical.

	    std::vector< double > vtemps(2);
	    std::vector< double > vdensities(2);
	    
	    vtemps[0] = temperature;
	    vtemps[1] = temperature;
	    vdensities[0] = density;
	    vdensities[1] = density;
	    
	    // Retrieve electron based heat capacities for each set of 
	    // (density, temperature) values.

	    std::vector< double > vCve(2);
	    vCve = spEospac->getElectronHeatCapacity( vtemps,
						      vdensities );

	    // Since the i=0 and i=1 tuples of density and temperature 
	    // are identical the two returned heat capacities should
	    // also match.

	    if ( match( vCve[0], vCve[1] ) )
		pass() << "getElectronHeatCapacity() test passed for "
		       << "vector state values.";
	    else
		fail() << "getElectronHeatCapacity() test failed for "
		       << "vector state values.";
	    
	    // This result should also match the scalar value
	    // calculated above.

	    heatCapacity =
		spEospac->getElectronHeatCapacity( 
		    temperature, density );

	    if ( match( vCve[0], heatCapacity ) )
		pass() << "getElectronHeatCapacity() test passed for "
		       << "vector state values.";
	    else
		fail() << "getElectronHeatCapacity() test failed for "
		       << "vector state values.";


	    // This feature has be removed from the public interface
	    // of cdi_eospac:

	    // The user can obtian information from EOSPAC by
	    // specifying the data type requested along with the
	    // (T,rho) tuple.  However, getdFdT() and getF() both
	    // expect vector argumetns for temperature and density AND 
	    // returns a vector of results.	    
	    
// 	    std::vector< double > altHeatCapacity = spEospac->getdFdT( 
// 		vtemps, vdensities, rtt_cdi_eospac::ES4enelc );

// 	    if ( match( altHeatCapacity, vCve ) )
// 		pass() << "getdFdT() test passed for "
// 		    "vector state values.";
// 	    else 
// 		fail() << "getdFdT() test failed for vector state values.";


	    // ---------------------- //
	    // Print the test result. //
	    // ---------------------- //
	    
	    if ( passed() ) {
		pass() << "All tests passed.";
		return "All tests passed.";
	    }
	    
	    return "Some tests failed.";
	    
	} // end of runTest()
    
    // ------------------------------------------ //
    // Compare Reference value to computed values //
    // ------------------------------------------ //
    
    bool tEospac::match( const double computedValue,
			 const double referenceValue ) const
	{
	    // Start by assuming that the two quantities match exactly.
	    bool em = true;
	    
	    // Compare items up to 10 digits of accuracy.
	    
	    const double TOL = 1.0e-10;
	    
	    // Calculate the absolute value of the relative difference between 
	    // the computed and reference values.
	    
// 	    std::cout.precision(16);
// 	    std::cout << "\t" << computedValue
// 		      << "\t" << referenceValue << std::endl;
	    
	    double reldiff = fabs( ( computedValue - referenceValue )
				   / referenceValue );
	    
	    // If the comparison fails then change the value of "em" return
	    // the result;
	    if ( reldiff > TOL )
		em = false;
	    
	    return em;    
	    
	} // end of tEospac::match( double, double )
    
    // ------------- //
    // Match vectors //
    // ------------- //
    
    bool tEospac::match( 
	const std::vector< double >& computedValue, 
	const std::vector< double >& referenceValue ) const
	{
	    // Start by assuming that the two quantities match exactly.
	    bool em = true;
	    
	    // Compare items up to 10 digits of accuracy.
	    const double TOL = 1.0e-10;
	    
	    // Test each item in the list
	    double reldiff = 0.0;
	    for ( int i=0; i<computedValue.size(); ++i )
		{
// 		    std::cout.precision(10);
// 		    std::cout << "\t" << computedValue[i]
// 			      << "\t" << referenceValue[i] << std::endl;
		    
		    reldiff = fabs( ( computedValue[i] - referenceValue[i] )
				    / referenceValue[i] );

		    // If the comparison fails then change the value of "em"
		    // and exit the loop.
		    
		    // DEBUG: must #include <iomanip>
		    //
		    // 	    std::cout << std::setprecision(14) << "   "
		    // 		      << computedValue[i] << "   "
		    // 		      << referenceValue[i] << "   "
		    // 		      << reldiff << std::endl;
		    
		    if ( reldiff > TOL )
			{
			    em = false;
			    break;
			}
		}
	    return em;
	} 
    
} // end namespace rtt_cdi_eospac_test

//---------------------------------------------------------------------------//
// end of tEospac.cc
//---------------------------------------------------------------------------//
