//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tDummyOpacity.cc
 * \author Kelly Thompson
 * \date   Fri Jan 5 15:59:00 2001
 * \brief  Implementation file for tDummyOpacity
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <string>
#include <vector>
#include <ostream>

#include "tDummyOpacity.hh"
#include "../Release.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

// Unit Test Frame Stuff
//----------------------------------------
namespace rtt_UnitTestFrame {
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
					   std::ostream& os_in ) {
	return rtt_dsxx::SP<TestApp> ( 
	    new rtt_cdi_dummy_opacity_test::tDummyOpacity( argc, argv, os_in ));
    }
} // end namespace rtt_UnitTestFrame

// tDummyOpacity Stuff
//--------------------------------------------------
namespace rtt_cdi_dummy_opacity_test {

using std::string;
using std::cout;
using std::endl;
using std::vector;
using rtt_dsxx::SP;

tDummyOpacity::tDummyOpacity( int argc, char *argv[], std::ostream& os_in )
    : rtt_UnitTestFrame::TestApp( argc, argv, os_in )
{
    os() << "Created tDummyOpacity" << endl;
}

string tDummyOpacity::version() const
{
    return rtt_cdi::release();
}

//===========================================================================
/*!
 * \brief This function runs a series of tests to exercise the
 *        functionality of the GandolfOpacity class.  These tests
 *        depend on the GandolfFile class but not on the Opacity
 *        class.
 *
 * \sa A GandolfFile object encapsulates all access to a single
 * IPCRESS file.  The user must specify the IPCRESS filename in the
 * constructor.  Many materials may be defined in a single IPCRESS
 * data file.  The GandolfOpacity object defines a single state of a
 * single IPCRESS material.  A GandolfOpacity object 
 * requires a GandolfFile object, an enumerator that specifies the
 * opacity model {rosseland or plank}, an enumerator that specifies
 * the reaction type {total, absorption or scattering} and an energy
 * policy model {Gray or Multigroup}.
 *
 */
//===========================================================================
string tDummyOpacity::runTest()
{
//     // Gandolf data filename (IPCRESS format required)
//     string op_data_file = "Al_BeCu.ipcress";
    
//     // ------------------------- //
//     // Create GandolfFile object //
//     // ------------------------- //
    
//     SP<rtt_cdi_gandolf::GandolfFile> spGFABC;

//     // Attempt to instantiate the object.
//     try
// 	{
// 	    spGFABC = new rtt_cdi_gandolf::GandolfFile( op_data_file ); 
// 	}
//     catch ( rtt_cdi_gandolf::GandolfException gerr )
// 	{
// 	    fail() << std::endl << "\t" << gerr.errorSummary();
// 	    return "Unable to create GanolfFile object. Test sequence aborted.";
// 	}

//     // If we make it here then spGFABC was successfully instantiated.
//     pass() << "SP to new GandolfFile object created for Al_BeCu.ipcress data.";

//     // Test the GandolfFile object.
//     if ( spGFABC->getDataFilename() == op_data_file )
// 	pass() << "GandolfFile object is now linked to the Al_BeCu.ipcress data file.";
//     else
// 	fail() << "GandolfFile object failed to link itself to the "
// 	       << "Al_BeCu.ipcress  data file.";

//     if ( spGFABC->getNumMaterials() == 2 )
// 	pass() << "The correct number of materials was found in the "
// 	       << "Al_BeCu.ipcress data file.";
//     else
// 	fail() << "spGFABC did not find the correct number of materials "
// 	       << "in the Al_BeCu.ipcress data file.";

//     // ---------------------- //
//     // Create Opacity object. //
//     // ---------------------- //

//     // Material identifier.  This data file has two materials: Al and
//     // BeCu.  Al has the id tag "10001".
//     const int matid=10001;
    
//     // SP< rtt_cdi_gandolf::GandolfOpacity<rtt_cdi_gandolf::Gray> > spOpacity_Al;
//     SP< rtt_cdi_gandolf::GrayOpacity > spOp_Al_rgt;
    
//     // Try to instantiate the Opacity object.
//     try
// 	{
// 	    spOp_Al_rgt = new rtt_cdi_gandolf::GrayOpacity( 
// 		spGFABC, 
// 		matid,
// 		rtt_cdi_gandolf::Rosseland, 
// 		rtt_cdi_gandolf::Total );
// 	}
//     catch ( rtt_cdi_gandolf::GandolfException gerr )
// 	// Alternatively, we could use:
// 	// catch ( rtt_cdi_gandolf::gkeysException gerr )
// 	// catch ( rtt_cdi_gandolf::gchgridsException gerr )
// 	// catch ( rtt_cdi_gandolf::ggetmgException gerr )
// 	// catch ( rtt_cdi_gandolf::ggetgrayException gerr )
// 	{
// 	    fail() << "Failed to create SP to new GandolfOpacity object for "
// 		   << "Al_BeCu.ipcress data."
// 		   << std::endl << "\t" << gerr.errorSummary();
// 	    return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
// 	}

//     // If we get here then the object was successfully instantiated.
//     pass() << "SP to new Opacity object created for Al_BeCu.ipcress data.";

    
//     // ----------------- //
//     // Gray Opacity Test //
//     // ----------------- //

//     double temperature = 0.1; // keV
//     double density = 27.0; // g/cm^3
//     double tabulatedGrayOpacity = 4271.7041147070677; // cm^2/g
    
//     if ( ! testGrayRosselandOpacityAccessorPassed( 
// 	spOp_Al_rgt, temperature, density, tabulatedGrayOpacity ) )
// 	return "testGrayRosselandOpacityAccessorPassed() failed.";
    
//     // --------------- //
//     // MG Opacity test //
//     // --------------- //

//     // Create a Multigroup Rosseland Total Opacity object (again for Al).
//     SP< rtt_cdi_gandolf::MultigroupOpacity > spOp_Al_rtmg;
    
//     // Try to instantiate the Opacity object.
//     try
// 	{
// 	    spOp_Al_rtmg = new rtt_cdi_gandolf::MultigroupOpacity( 
// 		spGFABC, 
// 		matid,
// 		rtt_cdi_gandolf::Rosseland, 
// 		rtt_cdi_gandolf::Total );
// 	}
//     catch ( rtt_cdi_gandolf::GandolfException gerr )
// 	{
// 	    fail() << "Failed to create SP to new GandolfOpacity object for "
// 		   << "Al_BeCu.ipcress data."
// 		   << std::endl << "\t" << gerr.errorSummary();
// 	    return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
// 	}

//   // Setup the test point.
//     temperature = 0.01; // keV
//     density = 2.0; // g/cm^3
    
//     // The solution to compare against:
//     int numGroups = 33;
//     double tabulatedMGOpacityArray[] =
//     {
// 	2.4935245299837247e+08,
// 	2.6666789027326573e+04,
// 	1.6270621515227660e+04,
// 	1.7634711671468522e+04,
// 	4.4999455359684442e+04,
// 	9.9917674335613032e+04,
// 	8.3261383385464113e+04,
// 	5.9742975304886764e+04,
// 	4.0373209722602740e+04,
// 	2.6156503146710609e+04,
// 	1.6356701105166874e+04,
// 	1.0007184686170869e+04,
// 	5.9763667878215247e+03,
// 	3.5203912630050986e+03,
// 	2.0765528559140448e+03,
// 	6.8550529299142445e+03,
// 	4.1257095227186965e+03,
// 	2.4199006949490426e+03,
// 	1.3894677080938793e+03,
// 	7.9046985091966621e+02,
// 	4.4088463936537232e+02,
// 	2.4514360684176387e+02,
// 	1.3541656611912146e+02,
// 	7.1828886317050177e+01,
// 	3.9793827527329107e+01,
// 	2.3312673181867030e+01,
// 	1.4879458895157605e+01,
// 	1.0862672679200283e+01,
// 	9.0590676798691288e+00,
// 	8.2841367649864175e+00,
// 	7.3809286930540363e+00,
// 	7.1057875403123791e+00,
// 	6.8907716134926735e+00
//     }; // KeV, numGroups entries.

//     vector<double> tabulatedMGOpacity(numGroups);
//     std::copy( tabulatedMGOpacityArray, 
// 	       tabulatedMGOpacityArray+numGroups,
// 	       tabulatedMGOpacity.begin() );    

//     if ( ! testMGRosselandOpacityAccessorPassed( 
// 	spOp_Al_rtmg, temperature, density, tabulatedMGOpacity ) )
// 	return "testGrayRosselandOpacityAccessorPassed() failed.";

//     // ----------------------------------------------- //
//     // Test the data file "analyticOpacities.ipcress"  //
//     // ----------------------------------------------- //

//      // -----------------------------------------------------------------
//      // The Opacities in this file are computed from the following
//      // analytic formula:
//      //     opacity = rho * T^4,
//      // rho is the density and T is the temperature.
//      //
//      // The grid in this data file has the following structure:
//      //    T   = { 0.1, 1.0, 10.0 } keV.
//      //    rho = { 0.1, 0.5, 1.0 } g/cm^3
//      //    E_bounds = { 0.01, 0.03, 0.07, 0.1, 0.3, 0.7, 1.0, 3.0, 7.0 
//      //                 10.0, 30.0, 70.0 100.0 } keV.
//      //-----------------------------------------------------------------

//      // Gandolf data filename (IPCRESS format required)
//      op_data_file = "analyticOpacities.ipcress";
     
//      // ------------------------- //
//      // Create GandolfFile object //
//      // ------------------------- //
     
//      // Create a smart pointer to a GandolfFile object
//      SP<rtt_cdi_gandolf::GandolfFile> spGFAnalytic;

//      // Try to instantiate the object.
//      try 
// 	 {
// 	     spGFAnalytic = new rtt_cdi_gandolf::GandolfFile( op_data_file ); 
// 	 }
//      catch ( rtt_cdi_gandolf::gmatidsException gerr)
// 	 {
// 	     fail() << std::endl << "\t" << gerr.errorSummary();
// 	     return "Unable to instantiate GandolfFile object spGFAnalytic.  Test sequence aborted.";
// 	 }
     
//      // If we make it here then spGFAnalytic was successfully instantiated.
//      pass() << "SP to new GandolfFile object created (spGFAnalytic).";
     
//      // Test the GandolfFile object.
//      if ( spGFAnalytic->getDataFilename() == op_data_file )
// 	 pass() << "GandolfFile object is now linked to the data file.";
//      else
// 	 fail() << "GandolfFile object failed to link itself to the data file.";

//      if ( spGFAnalytic->getNumMaterials() == 1 )
// 	 pass() << "The correct number of materials was found in the data file.";
//      else
// 	 fail() << "spGFAnalytic did not find the correct number of materials in the data file.";

//      // --------------------- //
//      // Create Opacity object //
//      // --------------------- //

//      // Create a smart pointer to an Opacity object.
//       SP<rtt_cdi_gandolf::GrayOpacity> spOp_Analytic_ragray;
     
//      // Try to instantiate the Opacity object.
//      try 
// 	 {
// 	     spOp_Analytic_ragray
// 		 = new rtt_cdi_gandolf::GrayOpacity( spGFAnalytic, 
// 						     matid,
// 						     rtt_cdi_gandolf::Rosseland,
// 						     rtt_cdi_gandolf::Absorption );
// 	 }
//       catch ( rtt_cdi_gandolf::GandolfException gerr )
// 	  // Alternatively, we could use:
// 	  // catch ( rtt_cdi_gandolf::gkeysException gerr )
// 	  // catch ( rtt_cdi_gandolf::gchgridsException gerr )
// 	  // catch ( rtt_cdi_gandolf::ggetmgException gerr )
// 	  // catch ( rtt_cdi_gandolf::ggetgrayException gerr )
// 	  {
// 	      fail() << "Failed to create SP to new GandolfOpacity object for "
// 		     << "Al_BeCu.ipcress data."
// 		     << std::endl << "\t" << gerr.errorSummary();
// 	      return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
// 	  }
      
//       // If we get here then the object was successfully instantiated.
//       pass() << "SP to new Opacity object created for analyticOpacities.ipcress.";

//       // ----------------- //
//       // Gray Opacity Test //
//       // ----------------- //
      
//       temperature = 10.0; // keV
//       density = 1.0; // g/cm^3
//       tabulatedGrayOpacity = density * pow( temperature, 4 );
      
//       if ( ! testGrayRosselandOpacityAccessorPassed( 
// 	  spOp_Analytic_ragray, temperature, density,
// 	  tabulatedGrayOpacity ) )
// 	  return "testGrayRosselandOpacityAccessorPassed() failed.";
     
//       //---------------- //
//       // MG Opacity test //
//       //---------------- //
      
//       typedef rtt_cdi_gandolf::MultigroupOpacity MgOp;
//       typedef rtt_cdi_gandolf::GrayOpacity GrOp;

//       // Create a smart pointer to an Opacity object.
//       SP<MgOp> spOp_Analytic_ramg;
     
//       // Try to instantiate the Opacity object.
//       try { spOp_Analytic_ramg
// 		= new MgOp( 
// 		    spGFAnalytic, matid, rtt_cdi_gandolf::Rosseland,
// 		    rtt_cdi_gandolf::Absorption );
//       } catch ( rtt_cdi_gandolf::GandolfException gerr ) {
// 	  fail() << "Failed to create SP to new GandolfOpacity object for "
// 		 << "Al_BeCu.ipcress data." << std::endl << "\t" << gerr.errorSummary();
// 	  return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
//       }
      
//       // If we get here then the object was successfully instantiated.
//       pass() << "SP to new Opacity object created for analyticOpacities.ipcress.";
      
//       // Set up the new test problem.
      
//       temperature = 0.3; // keV
//       density = 0.7; // g/cm^3
      
//       // This is the solution we compare against.
//       numGroups = 12;
//       tabulatedMGOpacity.resize(12);
//       for ( int i=0; i<numGroups; ++i )
// 	  tabulatedMGOpacity[i] = density * pow( temperature, 4 ); // cm^2/gm

//       if ( ! testMGRosselandOpacityAccessorPassed( 
// 	  spOp_Analytic_ramg, temperature, density, tabulatedMGOpacity ) )
// 	  return "testMGRosselandOpacityAccessorPassed() failed.";

//       // ------------------------------------------------------------ //
//       // Test the Plank routines using analyticOpacities.ipcress data //
//       // ------------------------------------------------------------ //
      
//       // The Opacities in this file are computed from the following
//       // analytic formula:
//       //     opacity = rho * T^4,
//       // rho is the density and T is the temperature.
      
//       // spGFAnalytic already points to the correct file so we don't repeat the 
//       // coding.
      
//       // Dito for spOpacityAnalytic.
      
//       // ----------------- //
//       // Gray Opacity Test //
//       // ----------------- //

//       // Create a smart pointer to an Opacity object.
//       SP<GrOp> spOp_Analytic_pgray;
     
//       // Try to instantiate the Opacity object.
//       try { spOp_Analytic_pgray
// 		= new GrOp( 
// 		    spGFAnalytic, matid, rtt_cdi_gandolf::Plank,
// 		    rtt_cdi_gandolf::Total );
//       } catch ( rtt_cdi_gandolf::GandolfException gerr ) {
// 	  fail() << "Failed to create SP to new GandolfOpacity object for "
// 		 << "Al_BeCu.ipcress data." << std::endl << "\t" << gerr.errorSummary();
// 	  return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
//       }
      
//       // If we get here then the object was successfully instantiated.
//       pass() << "SP to new Gray Plank Total Opacity object created for analyticOpacities.ipcress.";

//       // Setup the test problem.
      
//       temperature = 3.0; // keV
//       density = 0.7; // g/cm^3
//       double tabulatedValue = density * pow( temperature, 4 ); // cm^2/g
      
//       if ( ! testGrayPlankOpacityAccessorPassed( spOp_Analytic_pgray,
// 						 temperature, density,
// 						 tabulatedValue ) )
// 	  return "testGrayPlankOpacityAccessor() failed.";
      
//       // --------------- //
//       // MG Opacity test //
//       // --------------- //

//       // Create a smart pointer to an Opacity object.
//       SP<MgOp> spOp_Analytic_pmg;
     
//       // Try to instantiate the Opacity object.
//       try { spOp_Analytic_pmg
// 		= new MgOp( spGFAnalytic, matid,
// 			    rtt_cdi_gandolf::Plank,
// 			    rtt_cdi_gandolf::Total ); 
//       } catch ( rtt_cdi_gandolf::GandolfException gerr ) {
// 	  fail() << "Failed to create SP to new GandolfOpacity object for "
// 		 << "Al_BeCu.ipcress data." << std::endl << "\t" 
// 		 << gerr.errorSummary();
// 	  return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
//       }

//       // If we get here then the object was successfully instantiated.
//       pass() << "SP to new Multigroup Plank Total Opacity object created \n\t"
// 	     << "for \"analyticOpacities.ipcress.\"";

//       // Setup the test problem.

//       int ng=12;
//       tabulatedMGOpacity.resize( ng );
//       temperature = 0.4; // keV
//       density = 0.22; // g/cm^3
//       for ( int ig=0; ig<ng; ++ig )
// 	  tabulatedMGOpacity[ig] = density * pow( temperature, 4 ); // cm^2/g
      
//       // If this test fails then stop testing.
//       if ( ! testMGPlankOpacityAccessorPassed( spOp_Analytic_pmg, 
// 					       temperature, density,
// 					       tabulatedMGOpacity ) ) 
// 	  return "testMGPlankOpacityAccessor() failed.";
    
//     // ------------------------ //
//     // Access temperature grid. //
//     // ------------------------ //
    
//     testTemperatureGridAccessor( spOp_Analytic_pmg );
    
//     // ------------------------ //
//     // Access the density grid. //
//     // ------------------------ //
    
//     testDensityGridAccessor( spOp_Analytic_pmg );
    
//     // ----------------------------- //
//     // Access the energy boundaries. //
//     // ----------------------------- //
    
//     testEnergyBoundaryAccessor( spOp_Analytic_pmg );
    
//     // ------------------------------------------------------------ //
//     // Test alternate (vector-based) accessors for getGrayRosseland //
//     // ------------------------------------------------------------ //

//     // ---------------------- //
//     // Vector of temperatures //
//     // ---------------------- //

//     vector<double> vtemperature(2);
//     vtemperature[0] = 0.5; // keV
//     vtemperature[1] = 0.7; // keV
//     density = 0.35; // g/cm^3

//     vector<double> vtabulatedGrayOpacity( vtemperature.size() );
//     for ( int i=0; i< vtabulatedGrayOpacity.size(); ++i )
// 	vtabulatedGrayOpacity[i] = density * pow ( vtemperature[i], 4 );

//     if ( ! testGrayRosselandOpacityAccessorPassed( 
// 	spOp_Analytic_ragray, vtemperature, density,
// 	vtabulatedGrayOpacity ) )
// 	return "testGrayRosselandOpacityAccessorPassed() failed for a vector of temps.";

//     // ---------------------- //
//     // Vector of densities    //
//     // ---------------------- //
    
//     temperature = 0.3; //keV
//     vector<double> vdensity(3);
//     vdensity[0] = 0.2; // g/cm^3
//     vdensity[1] = 0.4; // g/cm^3
//     vdensity[2] = 0.6; // g/cm^3

//     vtabulatedGrayOpacity.resize( vdensity.size() );
//     for ( int i=0; i< vtabulatedGrayOpacity.size(); ++i )
// 	vtabulatedGrayOpacity[i] = vdensity[i] * pow ( temperature, 4 );

//     if ( ! testGrayRosselandOpacityAccessorPassed( 
// 	spOp_Analytic_ragray, temperature, vdensity,
// 	vtabulatedGrayOpacity ) )
// 	return "testGrayRosselandOpacityAccessorPassed() failed for a vector of densities.";
    
//     // ------------------------------------------------ //
//     // Vector of temperatures and a vector of densities //
//     // ------------------------------------------------ //

//     // This functionality was removed from GandolfOpacity because
//     // there is no way to ensure that the 2D container is filled in
//     // the correct order.  Consider using the STL-like accessor.

// //     int nt = vtemperature.size();
// //     int nd = vdensity.size();
// //     vtabulatedGrayOpacity.resize( nt*nd ); 
// //     for ( int it=0; it<nt; ++it )
// // 	for ( int id=0; id< nd; ++id )
// // 	    vtabulatedGrayOpacity[it*nd+id] = 
// // 		vdensity[id] * pow ( vtemperature[it], 4 );

// //     if ( ! testGrayRosselandOpacityAccessorPassed( 
// // 	spOp_Analytic_ragray, vtemperature, vdensity,
// // 	vtabulatedGrayOpacity ) )
// // 	return "testGrayRosselandOpacityAccessorPassed() failed for a vector of temperatures x densities.";

//     // -------------------------------------------------------- //
//     // Test alternate (vector-based) accessors for getGrayPlank //
//     // -------------------------------------------------------- //
    
//     // ---------------------- //
//     // Vector of temperatures //
//     // ---------------------- //
    
//     vtemperature.resize(2);
//     vtemperature[0] = 0.5; // keV
//     vtemperature[1] = 0.7; // keV
//     density = 0.35; // g/cm^3q
    
//     vtabulatedGrayOpacity.resize( vtemperature.size() );
//     for ( int i=0; i< vtabulatedGrayOpacity.size(); ++i )
// 	vtabulatedGrayOpacity[i] = density * pow ( vtemperature[i], 4 );

//     if ( ! testGrayPlankOpacityAccessorPassed( 
// 	spOp_Analytic_pgray, vtemperature, density,
// 	vtabulatedGrayOpacity ) )
// 	return "testGrayPlankOpacityAccessorPassed() failed for a vector of temps.";
   
//     // ------------------- //
//     // Vector of densities //
//     // ------------------- //

//     temperature = 0.3; //keV
//     vdensity.resize(3);
//     vdensity[0] = 0.2; // g/cm^3
//     vdensity[1] = 0.4; // g/cm^3
//     vdensity[2] = 0.6; // g/cm^3

//     vtabulatedGrayOpacity.resize( vdensity.size() );
//     for ( int i=0; i< vtabulatedGrayOpacity.size(); ++i )
// 	vtabulatedGrayOpacity[i] = vdensity[i] * pow ( temperature, 4 );

//     if ( ! testGrayPlankOpacityAccessorPassed( 
// 	spOp_Analytic_pgray, temperature, vdensity,
// 	vtabulatedGrayOpacity  ) )
// 	return "testGrayPlankOpacityAccessorPassed() failed for a vector of densities.";

//     // ------------------------------------------------ //
//     // Vector of temperatures and a vector of densities //
//     // ------------------------------------------------ //

//     // This functionality was removed from GandolfOpacity because
//     // there is no way to ensure that the 2D container is filled in
//     // the correct order. Consider using the STL-like accessor

// //     nt = vtemperature.size();
// //     nd = vdensity.size();
// //     vtabulatedGrayOpacity.resize( nt*nd ); 
// //     for ( int i=0; i< nt; ++i )
// // 	for ( int j=0; j< nd; ++j )
// // 	    vtabulatedGrayOpacity[i*nd+j] = 
// // 		vdensity[j] * pow ( vtemperature[i], 4 );
    
// //     if ( ! testGrayPlankOpacityAccessorPassed(
// // 	spOp_Analytic_pgray, vtemperature, vdensity,
// // 	vtabulatedGrayOpacity ) )
// //  	return "testGrayPlankOpacityAccessorPassed() failed for a vector of temperatures x densities.";

//     // ---------------------------------------------------------- //
//     // Test alternate (vector-based) accessors for getMGRosseland //
//     // ---------------------------------------------------------- //

//     // ---------------------- //
//     // Vector of temperatures //
//     // ---------------------- //

//      vtemperature.resize(2);
//      vtemperature[0] = 0.5; // keV
//      vtemperature[1] = 0.7; // keV
//      density = 0.35; // g/cm^3
//      ng = spOp_Analytic_ramg->getNumGroupBoundaries() - 1;

//      vector< vector< double > > vtabulatedMGOpacity( vtemperature.size() );
//      for ( int i=0; i< vtemperature.size(); ++i )
// 	 {
// 	     vtabulatedMGOpacity[i].resize(ng);
// 	     for ( int ig=0; ig<ng; ++ig )
// 		 vtabulatedMGOpacity[i][ig] = 
// 		     density * pow ( vtemperature[i], 4 );
// 	 }

//      if ( ! testMGRosselandOpacityAccessorPassed( spOp_Analytic_ramg, 
// 						  vtemperature, density,
// 						  vtabulatedMGOpacity ) )
// 	 return "testMGRosselandOpacityAccessorPassed() failed for a vector of temps.";

//      // ------------------- //
//      // Vector of densities //
//      // ------------------- //

//      vdensity.resize(2);
//      vdensity[0] = 0.3; // g/cm^3
//      vdensity[1] = 0.7; // g/cm^3
//      temperature = 7.0; // keV
//      ng = spOp_Analytic_ramg->getNumGroupBoundaries() - 1;

//      vtabulatedMGOpacity.resize( vdensity.size() );
//      for ( int i=0; i<vdensity.size(); ++i )
// 	 {
// 	     vtabulatedMGOpacity[i].resize(ng);
// 	     for ( int ig=0; ig<ng; ++ig )
// 		 vtabulatedMGOpacity[i][ig] = 
// 		     vdensity[i] * pow ( temperature, 4 );
// 	 }

//      if ( ! testMGRosselandOpacityAccessorPassed( 
// 	 spOp_Analytic_ramg, temperature, vdensity,
// 	 vtabulatedMGOpacity ) )
// 	 return "testMGRosselandOpacityAccessorPassed() failed for a vector of densities.";


//      // ------------------------------------------------ //
//      // Vector of temperatures and a vector of densities //
//      // ------------------------------------------------ //

//     // This functionality was removed from GandolfOpacity because
//     // there is no way to ensure that the 2D container is filled in
//     // the correct order.

// //      vtabulatedMGOpacity.resize( vtemperature.size() 
// // 				 * vdensity.size() );
// //      nd = vdensity.size();
// //      for ( int it=0; it<vtemperature.size(); ++it )
// // 	 for ( int id=0; id<nd; ++id )
// // 	     {
// // 		 vtabulatedMGOpacity[ it*nd + id ].resize(ng);
// // 		 for ( int ig=0; ig<ng; ++ig )
// // 		     vtabulatedMGOpacity[ it*nd + id ][ ig ] = 
// // 			 vdensity[id] * pow ( vtemperature[it], 4 );
// // 	     }

// //      if ( ! testMGRosselandOpacityAccessorPassed( 
// // 	 spOp_Analytic_ramg, vtemperature, vdensity,
// // 	 vtabulatedMGOpacity ) )
// // 	 return "testMGRosselandOpacityAccessorPassed() failed for a vector of temps and a vector of densities.";

//     // ------------------------------------------------------ //
//     // Test alternate (vector-based) accessors for getMGPlank //
//     // ------------------------------------------------------ //

//     // ---------------------- //
//     // Vector of temperatures //
//     // ---------------------- //

//      vtemperature.resize(2);
//      vtemperature[0] = 0.5; // keV
//      vtemperature[1] = 0.7; // keV
//      density = 0.35; // g/cm^3
//      ng = spOp_Analytic_pmg->getNumGroupBoundaries() - 1;

//      vtabulatedMGOpacity.resize( vtemperature.size() );
//      for ( int i=0; i< vtemperature.size(); ++i )
// 	 {
// 	     vtabulatedMGOpacity[i].resize(ng);
// 	     for ( int ig=0; ig<ng; ++ig )
// 		 vtabulatedMGOpacity[ i ][ ig ] = 
// 		     density * pow ( vtemperature[i], 4 );
// 	 }
     
//      if ( ! testMGPlankOpacityAccessorPassed( 
// 	 spOp_Analytic_pmg, vtemperature, density,
// 	 vtabulatedMGOpacity ) )
// 	 return "testMGPlankOpacityAccessorPassed() failed for a vector of temps.";

//      // ------------------- //
//      // Vector of densities //
//      // ------------------- //

//      vdensity.resize(2);
//      vdensity[0] = 0.3; // g/cm^3
//      vdensity[1] = 0.7; // g/cm^3
//      temperature = 7.0; // keV
//      ng = spOp_Analytic_pmg->getNumGroupBoundaries() - 1;

//      vtabulatedMGOpacity.resize( vdensity.size() );
//      for ( int i=0; i<vdensity.size(); ++i )
// 	 {
// 	     vtabulatedMGOpacity[i].resize(ng);
// 	     for ( int ig=0; ig<ng; ++ig )
// 		 vtabulatedMGOpacity[ i ][ ig ] = 
// 		     vdensity[i] * pow ( temperature, 4 );
// 	 }

//      if ( ! testMGPlankOpacityAccessorPassed( 
// 	 spOp_Analytic_pmg, temperature, vdensity,
// 	 vtabulatedMGOpacity ) )
// 	 return "testMGPlankOpacityAccessorPassed() failed for a vector of densities.";


//      // ------------------------------------------------ //
//      // Vector of temperatures and a vector of densities //
//      // ------------------------------------------------ //

//     // This functionality was removed from GandolfOpacity because
//     // there is no way to ensure that the 2D container is filled in
//     // the correct order.  Consider using the STL-like accessor.

// //      vtabulatedMGOpacity.resize( vtemperature.size() 
// // 				 * vdensity.size() );
// //      nd = vdensity.size();
// //      for ( int it=0; it<vtemperature.size(); ++it )
// // 	 for ( int id=0; id<nd; ++id )
// // 	     {
// // 		 vtabulatedMGOpacity[ it*nd + id ].resize(ng);
// // 		 for ( int ig=0; ig<ng; ++ig )
// // 		     vtabulatedMGOpacity[ it*nd + id][ ig ] = 
// // 			 vdensity[id] * pow ( vtemperature[it], 4 );
// // 	     }

// //      if ( ! testMGPlankOpacityAccessorPassed( 
// // 	 spOp_Analytic_pmg, vtemperature, vdensity,
// // 	 vtabulatedMGOpacity ) )
// // 	 return "testMGPlankOpacityAccessorPassed() failed for a vector of temps and a vector of densities.";

     
//      // ------------------------------------- //
//      // Test the STL-like getOpacity accessor //
//      // Using const iterators                 //
//      // ------------------------------------- //

//      vdensity.resize(6);
//      vtemperature.resize(6);

//      // (temperature,density) tuples.

//      vtemperature[0] = 0.5; // keV
//      vdensity[0] = 0.2; // g/cm^3

//      vtemperature[1] = 0.7; // keV
//      vdensity[1] = 0.2; // g/cm^3

//      vtemperature[2] = 0.5; // keV
//      vdensity[2] = 0.4; // g/cm^3

//      vtemperature[3] = 0.7; // keV
//      vdensity[3] = 0.4; // g/cm^3

//      vtemperature[4] = 0.5; // keV
//      vdensity[4] = 0.6; // g/cm^3

//      vtemperature[5] = 0.7; // keV
//      vdensity[5] = 0.6; // g/cm^3

//      const vector<double> cvdensity = vdensity;
//      const vector<double> cvtemperature = vtemperature;

//      int nt = cvtemperature.size();
//      int nd = cvdensity.size();
//      vtabulatedGrayOpacity.resize( nt ); 

//      for ( int i=0; i<nt; ++i )
// 	 vtabulatedGrayOpacity[i] = 
// 	     cvdensity[i] * pow ( cvtemperature[i], 4 );

//      vector<double> graOp(nt);

//      spOp_Analytic_ragray->getOpacity( cvtemperature.begin(),
// 				       cvtemperature.end(), 
// 				       cvdensity.begin(),
// 				       cvdensity.end(), 
// 				       graOp.begin() );
//      if ( match( graOp, vtabulatedGrayOpacity ) )
// 	 pass() << spOp_Analytic_ragray->getDataDescriptor()
// 		<< " opacity computation was good for "
// 		<< spOp_Analytic_ragray->getDataFilename()
// 		<< " (const-iterator accessor, temp x density).";
//      else
// 	 {
// 	     fail() << spOp_Analytic_ragray->getDataDescriptor()
// 		    << " opacity value is out of spec. for "
// 		    << spOp_Analytic_ragray->getDataFilename()
// 		    << " (non-const-iterator accessor, temp x density).";
// 	     return false;
// 	 }

//      // ------------------------------------- //
//      // Test the STL-like getOpacity accessor //
//      // Using non-const iterator              //
//      // ------------------------------------- //

//      spOp_Analytic_ragray->getOpacity( vtemperature.begin(),
// 				       vtemperature.end(), 
// 				       vdensity.begin(),
// 				       vdensity.end(), 
// 				       graOp.begin() );
//      if ( match( graOp, vtabulatedGrayOpacity ) )
// 	 pass() << spOp_Analytic_ragray->getDataDescriptor()
// 		<< " opacity computation was good for "
// 		<< spOp_Analytic_ragray->getDataFilename()
// 		<< " (non-const-iterator accessor, temp x density).";
//      else
// 	 {
// 	     fail() << spOp_Analytic_ragray->getDataDescriptor()
// 		    << " opacity value is out of spec. for "
// 		    << spOp_Analytic_ragray->getDataFilename()
// 		    << " (non-const-iterator accessor, temp x density).";
// 	     return false;
// 	 }

//      // ------------------------------------- //
//      // Test the STL-like getOpacity accessor //
//      // const iterator (temperature only)     //
//      // ------------------------------------- //

//      graOp.resize( nt );
//      vtabulatedGrayOpacity.resize( nt );
//      for ( int it=0; it<nt; ++it )
// 	 vtabulatedGrayOpacity[it] = density * pow( vtemperature[it], 4 );

//      spOp_Analytic_ragray->getOpacity( cvtemperature.begin(),
// 				       cvtemperature.end(), 
// 				       density,
// 				       graOp.begin() );
//      if ( match( graOp, vtabulatedGrayOpacity ) )
// 	 pass() << spOp_Analytic_ragray->getDataDescriptor()
// 		<< " opacity computation was good for "
// 		<< spOp_Analytic_ragray->getDataFilename()
// 		<< " (const iterator accessor, vtemps).";
//      else
// 	 {
// 	     fail() << spOp_Analytic_ragray->getDataDescriptor()
// 		    << " opacity value is out of spec. for "
// 		    << spOp_Analytic_ragray->getDataFilename()
// 		    << " (const iterator accessor, vtemps).";
// 	     return false;
// 	 }

//      // ------------------------------------- //
//      // Test the STL-like getOpacity accessor //
//      // const iterator ( density only)        //
//      // ------------------------------------- //

//      graOp.resize( nd );
//      vtabulatedGrayOpacity.resize( nd );
//      for ( int id=0; id<nd; ++id )
// 	 vtabulatedGrayOpacity[id] = vdensity[id] * pow( temperature, 4 );

//      spOp_Analytic_ragray->getOpacity( temperature,
// 				       cvdensity.begin(),
// 				       cvdensity.end(), 
// 				       graOp.begin() );
//      if ( match( graOp, vtabulatedGrayOpacity ) )
// 	 pass() << spOp_Analytic_ragray->getDataDescriptor()
// 		<< " opacity computation was good for "
// 		<< spOp_Analytic_ragray->getDataFilename()
// 		<< " (const iterator accessor, vdensity).";
//      else
// 	 {
// 	     fail() << spOp_Analytic_ragray->getDataDescriptor()
// 		    << " opacity value is out of spec. for "
// 		    << spOp_Analytic_ragray->getDataFilename()
// 		    << " (const iterator accessor, vdensity).";
// 	     return false;
// 	 }

















//      //     fail() << "Test STL-like accessors for multigroup data!!!";
















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

bool tDummyOpacity::match( const double computedValue,
			     const double referenceValue ) const
{
    // Start by assuming that the two quantities match exactly.
    bool em = true;

    // Compare items up to 10 digits of accuracy.
    const double TOL = 1.0e-10;

    // Calculate the absolute value of the relative difference between 
    // the computed and reference values.
    double reldiff = fabs( ( computedValue - referenceValue )
			   / referenceValue );
    
    // If the comparison fails then change the value of "em" return
    // the result;
    if ( reldiff > TOL )
	em = false;

    return em;    

} // end of tDummyOpacity::match( double, double )

bool tDummyOpacity::match( const vector<double>& computedValue, 
			     const vector<double>& referenceValue ) const
{
    // Start by assuming that the two quantities match exactly.
    bool em = true;

    // Compare items up to 10 digits of accuracy.
    const double TOL = 1.0e-10;

//     std::cout << "size = " << computedValue.size() << std::endl
// 	      << "size = " << referenceValue.size() << std::endl;

    // Test each item in the list
    double reldiff = 0.0;
    for ( int i=0; i<computedValue.size(); ++i )
	{
	    reldiff = fabs( ( computedValue[i] - referenceValue[i] )
			    / referenceValue[i] );
	    // If the comparison fails then change the value of "em"
	    // and exit the loop.
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
} // end of tDummyOpacity::match( vector<double>, vector<double> )

bool tDummyOpacity::match( 
    const vector< vector<double> >& computedValue, 
    const vector< vector<double> >& referenceValue ) const 
{
    // Start by assuming that the two quantities match exactly.
    bool em = true;

    // Compare items up to 10 digits of accuracy.
    const double TOL = 1.0e-10;

    // Test each item in the list
    double reldiff = 0.0;
    for ( int i=0; i<computedValue.size(); ++i )
	{
	    for ( int j=0; j<computedValue[i].size(); ++j )
		{
		    reldiff = fabs( ( computedValue[i][j] - referenceValue[i][j] )
				    / referenceValue[i][j] );
	    // If the comparison fails then change the value of "em"
	    // and exit the loop.
		    if ( reldiff > TOL ) {
			em = false; break; }
		}
	}
    return em;
} // end of tDummyOpacity::match( vector<double>, vector<double> )

// ---------------------------------------- //
// Test the gray Rosseland opacity accessor // 
// ---------------------------------------- //

template < class temperatureType, class densityType, class opType >
bool tDummyOpacity::testGrayRosselandOpacityAccessorPassed(
    const SP<rtt_dummy_opacity::GrayOpacity> spOpacity, 
    const temperatureType temperature, 
    const densityType density, 
    const opType tabulatedValue )
{
     opType grayRosselandOpacity
 	= spOpacity->getOpacity( temperature, density );

    // Make sure that the interpolated value matches previous
    // interpolations. 

    if ( match( grayRosselandOpacity, tabulatedValue ) )
	pass() << spOpacity->getDataDescriptor()
	       << " opacity computation was good for \n\t" 
	       << "\"" << spOpacity->getDataFilename() << "\" data."; 
    else
	{
	    fail() << spOpacity->getDataDescriptor()
		   << " opacity value is out of spec. for \n\t"
		   << "\"" << spOpacity->getDataFilename() << "\" data."; 
	    return false;
	}

    // If we get here then the test passed.
    return true;
}

// ---------------------------------------------- //
// Test the multigroup Rosseland opacity accessor //
// ---------------------------------------------- //

template < class temperatureType, class densityType, class TestValueType >
bool tDummyOpacity::testMGRosselandOpacityAccessorPassed(
    const SP<rtt_dummy_opacity::MultigroupOpacity> spOpacity,
    const temperatureType temperature, 
    const densityType density,
    const TestValueType tabulatedValues )
{
    // Interpolate the multigroup opacities.
    TestValueType mgRosselandOpacity( tabulatedValues.size() );
    mgRosselandOpacity
	= spOpacity->getOpacity( temperature, density );

    // Compare the interpolated value with previous interpolations:

    if ( match( mgRosselandOpacity, tabulatedValues ) )
 	pass() << spOpacity->getDataDescriptor()
	       << " opacity computation was\n\t"
	       << "good for the data obtained from the file " 
	       << "\"" << spOpacity->getDataFilename() << "\".";
    else 
	{
	    fail() << spOpacity->getDataDescriptor()
		   << " opacity computation failed for \n\t"
		   << "the data obtained from the file " 
		   << "\"" << spOpacity->getDataFilename() << "\".";
	    return false;
	}

    // If we get here the test has passed.
    return true;

} // end of tDummyOpacity::testMGRosselandOpacityAccessorPassed( )

// ------------------------------------ //
// Test the gray Plank opacity accessor //
// ------------------------------------ //
template < class tempT, class denT, class opT >
bool tDummyOpacity::testGrayPlankOpacityAccessorPassed(
    const SP<rtt_dummy_opacity::GrayOpacity> spOpacity,
    const tempT temperature, const denT density, const opT tabulatedValue )
{
    // Call the interpolation routine.
    opT grayPlankOpacity 
	= spOpacity->getOpacity( temperature, density );

    // Make sure that the interpolated value matches previous
    // interpolations. 
    if ( match( grayPlankOpacity, tabulatedValue ) )
	pass() << spOpacity->getDataDescriptor()
	       << " opacity computation was good for "
	       << spOpacity->getDataFilename();
    else 
	{
	    fail() << spOpacity->getDataDescriptor()
		   << "opacity value is out of spec. for "
		   << spOpacity->getDataFilename();
	    return false;
	}
    
    // If we get here then the test passed.
    return true;

} // end of bool tDummyOpacity::testGrayOpacityAccessorPassed( )

// ------------------------------------------ //
// Test the multigroup Plank Opacity accessor //
// ------------------------------------------ //

template < class temperatureType, class densityType, class TestValueType >
bool tDummyOpacity::testMGPlankOpacityAccessorPassed( 
    const SP<rtt_dummy_opacity::MultigroupOpacity> spOpacity,
    const temperatureType temperature,
    const densityType density,
    const TestValueType tabulatedValues )
{
    TestValueType mgPlankOpacity( tabulatedValues.size() );
    
    // Interpolate the multigroup opacities.
    try
	{
	    mgPlankOpacity
		= spOpacity->getOpacity( temperature, density );
	}
    catch ( rtt_dummy_opacity::ggetmgException gerr )
	{
	    fail() << std::endl << "\t" << gerr.errorSummary();
	    return false;
	}
    catch ( rtt_dummy_opacity::gkeysException gerr )
    	{
	    fail() << std::endl << "\t" << gerr.errorSummary();
	    return false;
	}
    catch ( rtt_dummy_opacity::GandolfException gerr )
	{
	    fail() << std::endl << "\t" << gerr.errorSummary();
	    return false;
	}

    // Compare the interpolated value with previous interpolations:
    if ( match( mgPlankOpacity, tabulatedValues ) )
	pass() << "Multigroup Plank Opacity computation was good\n\t"
	       << "for the data obtained from the file "
	       << "\"" << spOpacity->getDataFilename() << "\"";
    else
	fail() << "Multigroup Plank Opacity computation failed for\n\t"
	       << "the data obtained from the file " 
	       << "\"" << spOpacity->getDataFilename() << "\".";

    return true;

} // end of void tDummyOpacity::testMGOpacityAccessor( )

// ---------------------------------- //
// Test the Temperature Grid Accessor //
// ---------------------------------- //

void tDummyOpacity::testTemperatureGridAccessor( 
    const SP<rtt_dummy_opacity::MultigroupOpacity> spOpacity )
{
    // Read the temperature grid from the IPCRESS file.     
    vector<double> temps = spOpacity->getTemperatureGrid();
    
    // Verify that the size of the temperature grid looks right.  If
    // it is the right size then compare the temperature grid data to 
    // the data specified when we created the IPCRESS file using TOPS.
    if ( temps.size() == spOpacity->getNumTemperatures() &&
	 temps.size() == 3 )
	{
	    pass() << "The number of temperature points found in the data\n\t" 
		   << "grid matches the number returned by the\n\t"
		   << "getNumTemperatures() accessor.";
	    
	    // The grid specified by TOPS has 3 temperature points.
	    vector<double> temps_ref( temps.size() );
	    temps_ref[0] = 0.1;
	    temps_ref[1] = 1.0;
	    temps_ref[2] = 10.0;
	    
	    // Compare the grids.
	    if ( match( temps, temps_ref ) )
		pass() << "Temperature grid matches.";
	    else
		fail() << "Temperature grid did not match.";
	}
    else
	{
	    fail() << "The number of temperature points found in the data\n\t"
		   << "grid does not match the number returned by the\n\t"
		   << "getNumTemperatures() accessor.";
	    fail() << "Did not test the results returned by\n\t"
		   << "getTemperatureGrid().";
	}
    
} // end of tDummyOpacity::testTemperatureGridAccessor( )

// ------------------------------ //
// Test the Density Grid Accessor //
// ------------------------------ //

void tDummyOpacity::testDensityGridAccessor( 
    const SP<rtt_dummy_opacity::MultigroupOpacity> spOpacity )
{
    // Read the grid from the IPCRESS file.     
    vector<double> density = spOpacity->getDensityGrid();
    
    // Verify that the size of the density grid looks right.  If
    // it is the right size then compare the density grid data to 
    // the data specified when we created the IPCRESS file using TOPS.
    if ( density.size() == 3 &&
	 density.size() == spOpacity->getNumDensities() )
	{
	    pass() << "The number of density points found in the data\n\t"
		   << "grid matches the number returned by the\n\t"
		   << "getNumDensities() accessor.";
	    
	    // The grid specified by TOPS has 3 density points
	    vector<double> density_ref( density.size() );
	    density_ref[0] = 0.1;
	    density_ref[1] = 0.5;
	    density_ref[2] = 1.0;
	    
	    // Compare the grids.
	    if ( match( density, density_ref ) )
		pass() << "Density grid matches.";
	    else
		fail() << "Density grid did not match.";
	}
    else
	{
	    fail() << "The number of density points found in the data\n\t"
		   << "grid does not match the number returned by the\n\t"
		   << "getNumDensities() accessor.";
	    fail() << "Did not test the results returned by\n\t"  
		   << "getDensityGrid().";
	}
} // end of void tDummyOpacity::testDensityGridAccessor( )

// --------------------------------- //
// Test the Energy Boundary Accessor //
// --------------------------------- //

void tDummyOpacity::testEnergyBoundaryAccessor( 
    const SP<rtt_dummy_opacity::MultigroupOpacity> spOpacity )
{
    
     // Read the grid from the IPCRESS file.     
     vector<double> ebounds = spOpacity->getGroupBoundaries();

     // Verify that the size of the group boundary grid looks right.  If
     // it is the right size then compare the energy groups grid data to 
     // the data specified when we created the IPCRESS file using TOPS.
     if ( ebounds.size() == 13 &&
	  ebounds.size() == spOpacity->getNumGroupBoundaries() )
	 {
	     pass() << "The number of energy boundary points found in the data\n\t"
		    << "grid matches the number returned by the\n\t"
		    << "getNumGroupBoundaries() accessor.";

	     // The grid specified by TOPS has 13 energy boundaries.
	     vector<double> ebounds_ref(ebounds.size());
	     ebounds_ref[0] = 0.01;
	     ebounds_ref[1] = 0.03;
	     ebounds_ref[2] = 0.07;
	     ebounds_ref[3] = 0.1;
	     ebounds_ref[4] = 0.3;
	     ebounds_ref[5] = 0.7;
	     ebounds_ref[6] = 1.0;
	     ebounds_ref[7] = 3.0;
	     ebounds_ref[8] = 7.0;
	     ebounds_ref[9] = 10.0;
	     ebounds_ref[10] = 30.0;
	     ebounds_ref[11] = 70.0;
	     ebounds_ref[12] = 100.0;

	     // Compare the grids.
	     if ( match( ebounds, ebounds_ref ) )
		 pass() << "Energy group boundary grid matches.";
	     else
		 fail() << "Energy group boundary grid did not match.";
	     
	 }
     else
	 {
	     fail() << "The number of energy boundary points found in the data\n\t"
		    << "grid does not match the number returned by the\n\t"
		    << "getNumGroupBoundaries() accessor.";
	     fail() << "Did not test the results returned by\n\t"  
		    << "getGroupBoundaries().";
	 }
    
} // end of testEnergyBoundaryAccessor( SP<Opacity> spOpacity )

} // end namespace rtt_cdi_gandolf_test

//---------------------------------------------------------------------------//
//                            end of tDummyOpacity.cc
//---------------------------------------------------------------------------//
