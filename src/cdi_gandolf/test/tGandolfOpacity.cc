//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tGandolfOpacity.cc
 * \author Kelly Thompson
 * \date   Thu Jun 22 13:07:00 2000
 * \brief  Implementation file for tGandolfOpacity
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "tGandolfOpacity.hh"
#include "tGandolfOpacity.t.hh"
#include "../Release.hh"

#include "../GandolfGrayOpacity.hh"
#include "../GandolfGrayOpacity.t.hh"
#include "../GandolfMultigroupOpacity.hh"
#include "../GandolfMultigroupOpacity.t.hh"
#include "../GandolfFile.hh"
#include "../GandolfException.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

#include <vector>

// Unit Test Frame Stuff
//----------------------------------------
namespace rtt_UnitTestFrame {
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
					   std::ostream& os_in ) {
	return rtt_dsxx::SP<TestApp> ( 
	    new rtt_cdi_gandolf_test::tGandolfOpacity( argc, argv, os_in ));
    }
} // end namespace rtt_UnitTestFrame

// tGandolfOpacity Stuff
//--------------------------------------------------
namespace rtt_cdi_gandolf_test {

tGandolfOpacity::tGandolfOpacity( int argc, char *argv[], std::ostream& os_in )
    : rtt_UnitTestFrame::TestApp( argc, argv, os_in )
{
    os() << "Created tGandolfOpacity" << std::endl;
}

std::string tGandolfOpacity::version() const
{
    return rtt_cdi_gandolf::release();
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
std::string tGandolfOpacity::runTest()
{
    // Gandolf data filename (IPCRESS format required)
    std::string op_data_file = "Al_BeCu.ipcress";
    
    // ------------------------- //
    // Create GandolfFile object //
    // ------------------------- //
    
    rtt_dsxx::SP< rtt_cdi_gandolf::GandolfFile > spGFABC;

    // Attempt to instantiate the object.
    try
	{
	    spGFABC = new rtt_cdi_gandolf::GandolfFile( op_data_file ); 
	}
    catch ( const rtt_cdi_gandolf::GandolfException& GandError )
	{
	    fail() << std::endl << "\t" << GandError.errorSummary();
	    return "Unable to create GanolfFile object. Test sequence aborted.";
	}

    // If we make it here then spGFABC was successfully instantiated.
    pass() << "SP to new GandolfFile object created for Al_BeCu.ipcress data.";

    // Test the GandolfFile object.
    if ( spGFABC->getDataFilename() == op_data_file )
	pass() << "GandolfFile object is now linked to the Al_BeCu.ipcress data file.";
    else
	fail() << "GandolfFile object failed to link itself to the "
	       << "Al_BeCu.ipcress  data file.";

    if ( spGFABC->getNumMaterials() == 2 )
	pass() << "The correct number of materials was found in the "
	       << "Al_BeCu.ipcress data file.";
    else
	fail() << "spGFABC did not find the correct number of materials "
	       << "in the Al_BeCu.ipcress data file.";

    // ---------------------- //
    // Create Opacity object. //
    // ---------------------- //

    // Material identifier.  This data file has two materials: Al and
    // BeCu.  Al has the id tag "10001".
    const int matid=10001;

    // Try to instantiate the Opacity object. (Rosseland, Gray Total
    // for material 10001 in the IPCRESS file pointed to by spGFABC).
    rtt_dsxx::SP< rtt_cdi::GrayOpacity > spOp_Al_rgt;

    try
	{
	    spOp_Al_rgt = new rtt_cdi_gandolf::GandolfGrayOpacity( 
		spGFABC, matid,	rtt_cdi::Rosseland, rtt_cdi::Total );
	}
    catch ( const rtt_cdi_gandolf::GandolfException& GandError )
	// Alternatively, we could use:
	// catch ( rtt_cdi_gandolf::gkeysException GandError )
	// catch ( rtt_cdi_gandolf::gchgridsException GandError )
	// catch ( rtt_cdi_gandolf::ggetmgException GandError )
	// catch ( rtt_cdi_gandolf::ggetgrayException GandError )
	{
	    fail() << "Failed to create SP to new GandolfOpacity object for "
		   << "Al_BeCu.ipcress data."
		   << std::endl << "\t" << GandError.errorSummary();
	    return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
	}

    // If we get here then the object was successfully instantiated.
    pass() << "SP to new Opacity object created for Al_BeCu.ipcress data.";

    
    // ----------------- //
    // Gray Opacity Test //
    // ----------------- //

    double temperature = 0.1; // keV
    double density = 27.0; // g/cm^3
    double tabulatedGrayOpacity = 4271.7041147070677; // cm^2/g
    
    if ( ! opacityAccessorPassed( 
	spOp_Al_rgt, temperature, density, tabulatedGrayOpacity ) )
	return "opacityAccessorPassed() failed.";
    
    // --------------- //
    // MG Opacity test //
    // --------------- //

    // Create a Multigroup Rosseland Total Opacity object (again for Al).
    rtt_dsxx::SP< rtt_cdi::MultigroupOpacity > spOp_Al_rtmg;
    
    // Try to instantiate the Opacity object.
    try
	{
	    spOp_Al_rtmg = new rtt_cdi_gandolf::GandolfMultigroupOpacity( 
		spGFABC, 
		matid,
		rtt_cdi::Rosseland, 
		rtt_cdi::Total );
	}
    catch ( const rtt_cdi_gandolf::GandolfException& GandError )
	{
	    fail() << "Failed to create SP to new GandolfOpacity object for "
		   << "Al_BeCu.ipcress data."
		   << std::endl << "\t" << GandError.errorSummary();
	    return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
	}

  // Setup the test point.
    temperature = 0.01; // keV
    density = 2.0; // g/cm^3
    
    // The solution to compare against:
    int numGroups = 33;
    double tabulatedMGOpacityArray[] =
    {
	2.4935245299837247e+08,
	2.6666789027326573e+04,
	1.6270621515227660e+04,
	1.7634711671468522e+04,
	4.4999455359684442e+04,
	9.9917674335613032e+04,
	8.3261383385464113e+04,
	5.9742975304886764e+04,
	4.0373209722602740e+04,
	2.6156503146710609e+04,
	1.6356701105166874e+04,
	1.0007184686170869e+04,
	5.9763667878215247e+03,
	3.5203912630050986e+03,
	2.0765528559140448e+03,
	6.8550529299142445e+03,
	4.1257095227186965e+03,
	2.4199006949490426e+03,
	1.3894677080938793e+03,
	7.9046985091966621e+02,
	4.4088463936537232e+02,
	2.4514360684176387e+02,
	1.3541656611912146e+02,
	7.1828886317050177e+01,
	3.9793827527329107e+01,
	2.3312673181867030e+01,
	1.4879458895157605e+01,
	1.0862672679200283e+01,
	9.0590676798691288e+00,
	8.2841367649864175e+00,
	7.3809286930540363e+00,
	7.1057875403123791e+00,
	6.8907716134926735e+00
    }; // KeV, numGroups entries.

    std::vector< double > tabulatedMGOpacity( numGroups );
    std::copy( tabulatedMGOpacityArray, 
	       tabulatedMGOpacityArray+numGroups,
	       tabulatedMGOpacity.begin() );    

    if ( ! opacityAccessorPassed( 
	spOp_Al_rtmg, temperature, density, tabulatedMGOpacity ) )
	return "opacityAccessorPassed() failed.";

    // ----------------------------------------------- //
    // Test the data file "analyticOpacities.ipcress"  //
    // ----------------------------------------------- //

     // -----------------------------------------------------------------
     // The Opacities in this file are computed from the following
     // analytic formula:
     //     opacity = rho * T^4,
     // rho is the density and T is the temperature.
     //
     // The grid in this data file has the following structure:
     //    T   = { 0.1, 1.0, 10.0 } keV.
     //    rho = { 0.1, 0.5, 1.0 } g/cm^3
     //    E_bounds = { 0.01, 0.03, 0.07, 0.1, 0.3, 0.7, 1.0, 3.0, 7.0 
     //                 10.0, 30.0, 70.0 100.0 } keV.
     //-----------------------------------------------------------------

     // Gandolf data filename (IPCRESS format required)
     op_data_file = "analyticOpacities.ipcress";
     
     // ------------------------- //
     // Create GandolfFile object //
     // ------------------------- //
     
     // Create a smart pointer to a GandolfFile object
     rtt_dsxx::SP<rtt_cdi_gandolf::GandolfFile> spGFAnalytic;

     // Try to instantiate the object.
     try 
	 {
	     spGFAnalytic = new rtt_cdi_gandolf::GandolfFile( op_data_file ); 
	 }
     catch ( const rtt_cdi_gandolf::gmatidsException& GandError)
	 {
	     fail() << std::endl << "\t" << GandError.errorSummary();
	     return "Unable to instantiate GandolfFile object spGFAnalytic.  Test sequence aborted.";
	 }
     
     // If we make it here then spGFAnalytic was successfully instantiated.
     pass() << "SP to new GandolfFile object created (spGFAnalytic).";
     
     // Test the GandolfFile object.
     if ( spGFAnalytic->getDataFilename() == op_data_file )
	 pass() << "GandolfFile object is now linked to the data file.";
     else
	 fail() << "GandolfFile object failed to link itself to the data file.";

     if ( spGFAnalytic->getNumMaterials() == 1 )
	 pass() << "The correct number of materials was found in the data file.";
     else
	 fail() << "spGFAnalytic did not find the correct number of materials in the data file.";

     // --------------------- //
     // Create Opacity object //
     // --------------------- //

     // Create a smart pointer to an Opacity object.
     rtt_dsxx::SP< rtt_cdi::GrayOpacity > spOp_Analytic_ragray;
     
     // Try to instantiate the Opacity object.
     try 
	 {
	     spOp_Analytic_ragray
		 = new rtt_cdi_gandolf::GandolfGrayOpacity( 
		     spGFAnalytic, 
		     matid,
		     rtt_cdi::Rosseland,
		     rtt_cdi::Absorption );
	 }
      catch ( const rtt_cdi_gandolf::GandolfException& GandError )
	  // Alternatively, we could use:
	  // catch ( rtt_cdi_gandolf::gkeysException GandError )
	  // catch ( rtt_cdi_gandolf::gchgridsException GandError )
	  // catch ( rtt_cdi_gandolf::ggetmgException GandError )
	  // catch ( rtt_cdi_gandolf::ggetgrayException GandError )
	  {
	      fail() << "Failed to create SP to new GandolfOpacity object for "
		     << "Al_BeCu.ipcress data."
		     << std::endl << "\t" << GandError.errorSummary();
	      return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
	  }
      
      // If we get here then the object was successfully instantiated.
      pass() << "SP to new Opacity object created for analyticOpacities.ipcress.";

      // ----------------- //
      // Gray Opacity Test //
      // ----------------- //
      
      temperature = 10.0; // keV
      density = 1.0; // g/cm^3
      tabulatedGrayOpacity = density * pow( temperature, 4 );
      
      if ( ! opacityAccessorPassed( 
	  spOp_Analytic_ragray, temperature, density,
	  tabulatedGrayOpacity ) )
	  return "opacityAccessorPassed() failed.";
     
      //---------------- //
      // MG Opacity test //
      //---------------- //
      
      // Create a smart pointer to an Opacity object.
      rtt_dsxx::SP< rtt_cdi::MultigroupOpacity > spOp_Analytic_ramg;
     
      // Try to instantiate the Opacity object.
      try { spOp_Analytic_ramg
		= new rtt_cdi_gandolf::GandolfMultigroupOpacity( 
		    spGFAnalytic, 
		    matid, 
		    rtt_cdi::Rosseland,
		    rtt_cdi::Absorption );
      } catch ( const rtt_cdi_gandolf::GandolfException& GandError ) {
	  fail() << "Failed to create SP to new GandolfOpacity object for "
		 << "Al_BeCu.ipcress data." << std::endl << "\t" << GandError.errorSummary();
	  return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
      }
      
      // If we get here then the object was successfully instantiated.
      pass() << "SP to new Opacity object created for analyticOpacities.ipcress.";
      
      // Set up the new test problem.
      
      temperature = 0.3; // keV
      density = 0.7; // g/cm^3
      
      // This is the solution we compare against.
      numGroups = 12;
      tabulatedMGOpacity.resize(12);
      for ( int i=0; i<numGroups; ++i )
	  tabulatedMGOpacity[i] = density * pow( temperature, 4 ); // cm^2/gm
      
      if ( ! opacityAccessorPassed( 
	  spOp_Analytic_ramg, temperature, density, tabulatedMGOpacity ) )
	  return "opacityAccessorPassed() failed.";
      
      // ------------------------------------------------------------ //
      // Test the Plank routines using analyticOpacities.ipcress data //
      // ------------------------------------------------------------ //
      
      // The Opacities in this file are computed from the following
      // analytic formula:
      //     opacity = rho * T^4,
      // rho is the density and T is the temperature.
      
      // spGFAnalytic already points to the correct file so we don't repeat the 
      // coding.
      
      // Dito for spOpacityAnalytic.
      
      // ----------------- //
      // Gray Opacity Test //
      // ----------------- //

      // Create a smart pointer to an Opacity object.
      rtt_dsxx::SP< rtt_cdi::GrayOpacity > spOp_Analytic_pgray;
     
      // Try to instantiate the Opacity object.
      try { spOp_Analytic_pgray
		= new rtt_cdi_gandolf::GandolfGrayOpacity( 
		    spGFAnalytic, 
		    matid, 
		    rtt_cdi::Plank,
		    rtt_cdi::Total );
      } catch ( const rtt_cdi_gandolf::GandolfException& GandError ) {
	  fail() << "Failed to create SP to new GandolfOpacity object for "
		 << "Al_BeCu.ipcress data." << std::endl << "\t" << GandError.errorSummary();
	  return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
      }
      
      // If we get here then the object was successfully instantiated.
      pass() << "SP to new Gray Plank Total Opacity object created for analyticOpacities.ipcress.";

      // Setup the test problem.
      
      temperature = 3.0; // keV
      density = 0.7; // g/cm^3
      double tabulatedValue = density * pow( temperature, 4 ); // cm^2/g
      
      if ( ! opacityAccessorPassed( spOp_Analytic_pgray,
						 temperature, density,
						 tabulatedValue ) )
	  return "testGrayPlankOpacityAccessor() failed.";
      
      // --------------- //
      // MG Opacity test //
      // --------------- //

      // Create a smart pointer to an Opacity object.
      rtt_dsxx::SP< rtt_cdi::MultigroupOpacity > spOp_Analytic_pmg;
     
      // Try to instantiate the Opacity object.
      try { spOp_Analytic_pmg
		= new rtt_cdi_gandolf::GandolfMultigroupOpacity( 
		    spGFAnalytic, 
		    matid,
		    rtt_cdi::Plank,
		    rtt_cdi::Total ); 
      } catch ( const rtt_cdi_gandolf::GandolfException& GandError ) {
	  fail() << "Failed to create SP to new GandolfOpacity object for "
		 << "Al_BeCu.ipcress data." << std::endl << "\t" 
		 << GandError.errorSummary();
	  return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
      }

      // If we get here then the object was successfully instantiated.
      pass() << "SP to new Multigroup Plank Total Opacity object created \n\t"
	     << "for \"analyticOpacities.ipcress.\"";

      // Setup the test problem.
      
      int ng=12;
      tabulatedMGOpacity.resize( ng );
      temperature = 0.4; // keV
      density = 0.22; // g/cm^3
      for ( int ig=0; ig<ng; ++ig )
	  tabulatedMGOpacity[ig] = density * pow( temperature, 4 ); // cm^2/g
      
      // If this test fails then stop testing.
      if ( ! opacityAccessorPassed( spOp_Analytic_pmg, 
				    temperature, density,
				    tabulatedMGOpacity ) ) 
	  return "testMGPlankOpacityAccessor() failed.";
    
      // ------------------------ //
      // Access temperature grid. //
      // ------------------------ //
      
      testTemperatureGridAccessor( spOp_Analytic_pmg );
      
      // ------------------------ //
      // Access the density grid. //
      // ------------------------ //
      
      testDensityGridAccessor( spOp_Analytic_pmg );
    
      // ----------------------------- //
      // Access the energy boundaries. //
      // ----------------------------- //
      
      testEnergyBoundaryAccessor( spOp_Analytic_pmg );
      
      // ------------------------------------------------------------ //
      // Test alternate (vector-based) accessors for getGrayRosseland //
      // ------------------------------------------------------------ //
      
      // ---------------------- //
      // Vector of temperatures //
      // ---------------------- //
      
      std::vector<double> vtemperature(2);
      vtemperature[0] = 0.5; // keV
      vtemperature[1] = 0.7; // keV
      density = 0.35; // g/cm^3
      
      std::vector<double> vtabulatedGrayOpacity( vtemperature.size() );
      for ( int i=0; i< vtabulatedGrayOpacity.size(); ++i )
	  vtabulatedGrayOpacity[i] = density * pow ( vtemperature[i], 4 );
      
      if ( ! opacityAccessorPassed( 
	  spOp_Analytic_ragray, vtemperature, density,
	  vtabulatedGrayOpacity ) )
	  return "opacityAccessorPassed() failed for a vector of temps.";

    // ---------------------- //
    // Vector of densities    //
    // ---------------------- //
    
    temperature = 0.3; //keV
    std::vector<double> vdensity(3);
    vdensity[0] = 0.2; // g/cm^3
    vdensity[1] = 0.4; // g/cm^3
    vdensity[2] = 0.6; // g/cm^3

    vtabulatedGrayOpacity.resize( vdensity.size() );
    for ( int i=0; i< vtabulatedGrayOpacity.size(); ++i )
	vtabulatedGrayOpacity[i] = vdensity[i] * pow ( temperature, 4 );

    if ( ! opacityAccessorPassed( 
	spOp_Analytic_ragray, temperature, vdensity,
	vtabulatedGrayOpacity ) )
	return "opacityAccessorPassed() failed for a vector of densities.";
    
    // -------------------------------------------------------- //
    // Test alternate (vector-based) accessors for getGrayPlank //
    // -------------------------------------------------------- //
    
    // ---------------------- //
    // Vector of temperatures //
    // ---------------------- //
    
    vtemperature.resize(2);
    vtemperature[0] = 0.5; // keV
    vtemperature[1] = 0.7; // keV
    density = 0.35; // g/cm^3q
    
    vtabulatedGrayOpacity.resize( vtemperature.size() );
    for ( int i=0; i< vtabulatedGrayOpacity.size(); ++i )
	vtabulatedGrayOpacity[i] = density * pow ( vtemperature[i], 4 );

    if ( ! opacityAccessorPassed( 
	spOp_Analytic_pgray, vtemperature, density,
	vtabulatedGrayOpacity ) )
	return "opacityAccessorPassed() failed for a vector of temps.";
   
    // ------------------- //
    // Vector of densities //
    // ------------------- //

    temperature = 0.3; //keV
    vdensity.resize(3);
    vdensity[0] = 0.2; // g/cm^3
    vdensity[1] = 0.4; // g/cm^3
    vdensity[2] = 0.6; // g/cm^3

    vtabulatedGrayOpacity.resize( vdensity.size() );
    for ( int i=0; i< vtabulatedGrayOpacity.size(); ++i )
	vtabulatedGrayOpacity[i] = vdensity[i] * pow ( temperature, 4 );

    if ( ! opacityAccessorPassed( 
	spOp_Analytic_pgray, temperature, vdensity,
	vtabulatedGrayOpacity  ) )
	return "opacityAccessorPassed() failed for a vector of densities.";

    // ---------------------------------------------------------- //
    // Test alternate (vector-based) accessors for getMGRosseland //
    // ---------------------------------------------------------- //

    // ---------------------- //
    // Vector of temperatures //
    // ---------------------- //

     vtemperature.resize(2);
     vtemperature[0] = 0.5; // keV
     vtemperature[1] = 0.7; // keV
     density = 0.35; // g/cm^3
     ng = spOp_Analytic_ramg->getNumGroupBoundaries() - 1;

     std::vector< std::vector< double > > vtabulatedMGOpacity( vtemperature.size() );
     for ( int i=0; i< vtemperature.size(); ++i )
	 {
	     vtabulatedMGOpacity[i].resize(ng);
	     for ( int ig=0; ig<ng; ++ig )
		 vtabulatedMGOpacity[i][ig] = 
		     density * pow ( vtemperature[i], 4 );
	 }

     if ( ! opacityAccessorPassed( spOp_Analytic_ramg, 
				   vtemperature, density,
				   vtabulatedMGOpacity ) )
	 return "opacityAccessorPassed() failed for a vector of temps.";

     // ------------------- //
     // Vector of densities //
     // ------------------- //

     vdensity.resize(2);
     vdensity[0] = 0.3; // g/cm^3
     vdensity[1] = 0.7; // g/cm^3
     temperature = 7.0; // keV
     ng = spOp_Analytic_ramg->getNumGroupBoundaries() - 1;

     vtabulatedMGOpacity.resize( vdensity.size() );
     for ( int i=0; i<vdensity.size(); ++i )
	 {
	     vtabulatedMGOpacity[i].resize(ng);
	     for ( int ig=0; ig<ng; ++ig )
		 vtabulatedMGOpacity[i][ig] = 
		     vdensity[i] * pow ( temperature, 4 );
	 }

     if ( ! opacityAccessorPassed( 
	 spOp_Analytic_ramg, temperature, vdensity,
	 vtabulatedMGOpacity ) )
	 return "opacityAccessorPassed() failed for a vector of densities.";

    // ------------------------------------------------------ //
    // Test alternate (vector-based) accessors for getMGPlank //
    // ------------------------------------------------------ //

    // ---------------------- //
    // Vector of temperatures //
    // ---------------------- //

     vtemperature.resize(2);
     vtemperature[0] = 0.5; // keV
     vtemperature[1] = 0.7; // keV
     density = 0.35; // g/cm^3
     ng = spOp_Analytic_pmg->getNumGroupBoundaries() - 1;

     vtabulatedMGOpacity.resize( vtemperature.size() );
     for ( int i=0; i< vtemperature.size(); ++i )
	 {
	     vtabulatedMGOpacity[i].resize(ng);
	     for ( int ig=0; ig<ng; ++ig )
		 vtabulatedMGOpacity[ i ][ ig ] = 
		     density * pow ( vtemperature[i], 4 );
	 }
     
     if ( ! opacityAccessorPassed( 
	 spOp_Analytic_pmg, vtemperature, density,
	 vtabulatedMGOpacity ) )
	 return "opacityAccessorPassed() failed for a vector of temps.";

     // ------------------- //
     // Vector of densities //
     // ------------------- //

     vdensity.resize(2);
     vdensity[0] = 0.3; // g/cm^3
     vdensity[1] = 0.7; // g/cm^3
     temperature = 7.0; // keV
     ng = spOp_Analytic_pmg->getNumGroupBoundaries() - 1;

     vtabulatedMGOpacity.resize( vdensity.size() );
     for ( int i=0; i<vdensity.size(); ++i )
	 {
	     vtabulatedMGOpacity[i].resize(ng);
	     for ( int ig=0; ig<ng; ++ig )
		 vtabulatedMGOpacity[ i ][ ig ] = 
		     vdensity[i] * pow ( temperature, 4 );
	 }

     if ( ! opacityAccessorPassed( 
	 spOp_Analytic_pmg, temperature, vdensity,
	 vtabulatedMGOpacity ) )
	 return "opacityAccessorPassed() failed for a vector of densities.";

     // -------------------------------------- //
     // Test the STL-like getOpacity accessor  //
     // Using const iterators for Gray objects //
     // -------------------------------------- //

     // These accessors are only available in GandolfOpacity objects
     // so the SP must be templated on GandolfGrayOpacity and not on
     // cdi/GrayOpacity. 
     
     // Create a new smart pointer to a GandolfGrayOpacity object.
     rtt_dsxx::SP< rtt_cdi_gandolf::GandolfGrayOpacity >
	 spGGOp_Analytic_ra;
     
     // try to instantiate the Opacity object.
     try
	 {
	     spGGOp_Analytic_ra = new
		 rtt_cdi_gandolf::GandolfGrayOpacity(
		     spGFAnalytic,
		     matid,
		     rtt_cdi::Rosseland,
		     rtt_cdi::Absorption );
	 }
     catch ( const rtt_cdi_gandolf::GandolfException& GandError )
	 {
	     fail() << "Failed to create SP to new GandolfGrayOpacity object for "
		    << "\n\tAl_BeCu.ipcress data (SP not templated on cdi/GrayOpacity)."
		    << std::endl << "\t" << GandError.errorSummary();
	     return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
	 }
     
     // If we get here then the object was successfully instantiated.
     pass() << "SP to new Opacity object created for analyticOpacities.ipcress.";
     
     // Setup the temperature and density parameters for this test.
     vdensity.resize(6);
     vtemperature.resize(6);
     
     // (temperature,density) tuples.
     
     vtemperature[0] = 0.5; // keV
     vdensity[0] = 0.2; // g/cm^3

     vtemperature[1] = 0.7; // keV
     vdensity[1] = 0.2; // g/cm^3

     vtemperature[2] = 0.5; // keV
     vdensity[2] = 0.4; // g/cm^3

     vtemperature[3] = 0.7; // keV
     vdensity[3] = 0.4; // g/cm^3

     vtemperature[4] = 0.5; // keV
     vdensity[4] = 0.6; // g/cm^3

     vtemperature[5] = 0.7; // keV
     vdensity[5] = 0.6; // g/cm^3

     // we want to test the const_iterator version of getOpacity() so
     // we need to create const vectors with the tuple data.
     const std::vector<double> cvdensity = vdensity;
     const std::vector<double> cvtemperature = vtemperature;

     int nt = cvtemperature.size();
     int nd = cvdensity.size();

     // Here is the reference solution
     vtabulatedGrayOpacity.resize( nt ); 
     for ( int i=0; i<nt; ++i )
	 vtabulatedGrayOpacity[i] = 
	     cvdensity[i] * pow ( cvtemperature[i], 4 );

     // Here is the solution from Gandolf
     std::vector< double > graOp(nt);
     spGGOp_Analytic_ra->getOpacity( cvtemperature.begin(),
				     cvtemperature.end(), 
				     cvdensity.begin(),
				     cvdensity.end(), 
				     graOp.begin() );
     if ( match( graOp, vtabulatedGrayOpacity ) )
	 pass() << spGGOp_Analytic_ra->getDataDescriptor()
		<< " opacity computation was good for \n\t"
		<< spGGOp_Analytic_ra->getDataFilename()
		<< " (const-iterator accessor, temp x density).";
     else
	 {
	     fail() << spGGOp_Analytic_ra->getDataDescriptor()
		    << " opacity value is out of spec. for \n\t"
		    << spGGOp_Analytic_ra->getDataFilename()
		    << " (non-const-iterator accessor, temp x density).";
	     return false;
	 }

     // ------------------------------------- //
     // Test the STL-like getOpacity accessor //
     // Using non-const iterator              //
     // ------------------------------------- //

     spGGOp_Analytic_ra->getOpacity( vtemperature.begin(),
				     vtemperature.end(), 
				     vdensity.begin(),
				     vdensity.end(), 
				     graOp.begin() );
     if ( match( graOp, vtabulatedGrayOpacity ) )
	 pass() << spGGOp_Analytic_ra->getDataDescriptor()
		<< " opacity computation was good for \n\t"
		<< spGGOp_Analytic_ra->getDataFilename()
		<< " (non-const-iterator accessor, temp x density).";
     else
	 {
	     fail() << spGGOp_Analytic_ra->getDataDescriptor()
		    << " opacity value is out of spec. for \n\t"
		    << spGGOp_Analytic_ra->getDataFilename()
		    << " (non-const-iterator accessor, temp x density).";
	     return false;
	 }
     
     // ------------------------------------- //
     // Test the STL-like getOpacity accessor //
     // const iterator (temperature only)     //
     // ------------------------------------- //
     
     graOp.resize( nt );
     vtabulatedGrayOpacity.resize( nt );
     for ( int it=0; it<nt; ++it )
	 vtabulatedGrayOpacity[it] = density * pow( vtemperature[it], 4 );

     spGGOp_Analytic_ra->getOpacity( cvtemperature.begin(),
				     cvtemperature.end(), 
				     density,
				     graOp.begin() );
     if ( match( graOp, vtabulatedGrayOpacity ) )
	 pass() << spGGOp_Analytic_ra->getDataDescriptor()
		<< " opacity computation was good for \n\t"
		<< spGGOp_Analytic_ra->getDataFilename()
		<< " (const iterator accessor, vtemps).";
     else
	 {
	     fail() << spGGOp_Analytic_ra->getDataDescriptor()
		    << " opacity value is out of spec. for \n\t"
		    << spGGOp_Analytic_ra->getDataFilename()
		    << " (const iterator accessor, vtemps).";
	     return false;
	 }

     // ------------------------------------- //
     // Test the STL-like getOpacity accessor //
     // const iterator ( density only)        //
     // ------------------------------------- //

     graOp.resize( nd );
     vtabulatedGrayOpacity.resize( nd );
     for ( int id=0; id<nd; ++id )
	 vtabulatedGrayOpacity[id] = vdensity[id] * pow( temperature, 4 );
     
     spGGOp_Analytic_ra->getOpacity( temperature,
				     cvdensity.begin(),
				     cvdensity.end(), 
				     graOp.begin() );
     if ( match( graOp, vtabulatedGrayOpacity ) )
	 pass() << spGGOp_Analytic_ra->getDataDescriptor()
		<< " opacity computation was good for \n\t"
		<< spGGOp_Analytic_ra->getDataFilename()
		<< " (const iterator accessor, vdensity).";
     else
	 {
	     fail() << spGGOp_Analytic_ra->getDataDescriptor()
		    << " opacity value is out of spec. for \n\t"
		    << spGGOp_Analytic_ra->getDataFilename()
		    << " (const iterator accessor, vdensity).";
	     return false;
	 }


     // -------------------------------------- //
     // Test the STL-like getOpacity accessor  //
     // Using const iterators for MG objects   //
     // -------------------------------------- //

     // These accessors are only available in GandolfOpacity objects
     // so the SP must be templated on GandolfMultigroupOpacity and not on
     // cdi/MultigroupOpacity. 
     
     // Create a new smart pointer to a GandolfGrayOpacity object.
     rtt_dsxx::SP< rtt_cdi_gandolf::GandolfMultigroupOpacity >
	 spGMGOp_Analytic_ra;
     
     // try to instantiate the Opacity object.
     try
	 {
	     spGMGOp_Analytic_ra = new
		 rtt_cdi_gandolf::GandolfMultigroupOpacity(
		     spGFAnalytic,
		     matid,
		     rtt_cdi::Rosseland,
		     rtt_cdi::Absorption );
	 }
     catch ( const rtt_cdi_gandolf::GandolfException& GandError )
	 {
	     fail() << "Failed to create SP to new GandolfGrayOpacity object for \n\t"
		    << "Al_BeCu.ipcress data (SP not templated on cdi/GrayOpacity)."
		    << std::endl << "\t" << GandError.errorSummary();
	     return "Unable to instantiate GandolfOpacity object.  Test sequence aborted.";
	 }

     // If we get here then the object was successfully instantiated.
     pass() << "SP to new Opacity object created for analyticOpacities.ipcress.";

     // Here is the reference solution
     ng = spGMGOp_Analytic_ra->getNumGroupBoundaries() - 1;
     std::vector< double > vtabulatedOpacity( ng * nt );
     
     for ( int i=0; i<nt; ++i )
	 for ( int ig=0; ig<ng; ++ig )
	     vtabulatedOpacity[ i*ng + ig ] = 
		 cvdensity[i] * pow ( cvtemperature[i], 4 );

     // Here is the solution from Gandolf
     std::vector< double > mgOp( nt*ng );
     spGMGOp_Analytic_ra->getOpacity( cvtemperature.begin(),
				      cvtemperature.end(), 
				      cvdensity.begin(),
				      cvdensity.end(), 
				      mgOp.begin() );

     if ( match( mgOp, vtabulatedOpacity ) )
	 pass() << spGMGOp_Analytic_ra->getDataDescriptor()
		<< " opacity computation was good for \n\t"
		<< spGMGOp_Analytic_ra->getDataFilename()
		<< " (const-iterator accessor, temp x density).";
     else
	 {
	     fail() << spGMGOp_Analytic_ra->getDataDescriptor()
		    << " opacity value is out of spec. for \n\t"
		    << spGMGOp_Analytic_ra->getDataFilename()
		    << " (non-const-iterator accessor, temp x density).";
	     return false;
	 }

     // --------------------------------------- //
     // Test the STL-like getOpacity accessor   //
     // Using non-const iterator for MG objects //
     // --------------------------------------- //

     // clear old data
     for ( int i=0; i<nt*ng; ++i) mgOp[i]=0.0;
     
     // use Gandolf to obtain new data
     spGMGOp_Analytic_ra->getOpacity( vtemperature.begin(),
				      vtemperature.end(), 
				      vdensity.begin(),
				      vdensity.end(), 
				      mgOp.begin() );

     // compare the results to the reference solution and report our
     // findings. 
     if ( match( mgOp, vtabulatedOpacity ) )
	 pass() << spGMGOp_Analytic_ra->getDataDescriptor()
		<< " opacity computation was good for \n\t"
		<< spGMGOp_Analytic_ra->getDataFilename()
		<< " (non-const-iterator accessor, temp x density).";
     else
	 {
	     fail() << spGMGOp_Analytic_ra->getDataDescriptor()
		    << " opacity value is out of spec. for \n\t"
		    << spGMGOp_Analytic_ra->getDataFilename()
		    << " (non-const-iterator accessor, temp x density).";
	     return false;
	 }

     // ------------------------------------------------ //
     // Test the STL-like getOpacity accessor            //
     // const iterator (temperature only) for MG data    //
     // ------------------------------------------------ //
     
     // clear old data
     for ( int i=0; i<nt*ng; ++i) mgOp[i]=0.0;

     // Calculate the reference solution.
     for ( int it=0; it<nt; ++it )
	 for ( int ig=0; ig<ng; ++ig )
	     vtabulatedOpacity[ it*ng + ig ] 
		 = density * pow( vtemperature[it], 4 );

     // Obtain new solution
     spGMGOp_Analytic_ra->getOpacity( cvtemperature.begin(),
				      cvtemperature.end(), 
				      density,
				      mgOp.begin() );

     // Compare solutions and report the results.
     if ( match( mgOp, vtabulatedOpacity ) )
	 pass() << spGGOp_Analytic_ra->getDataDescriptor()
		<< " opacity computation was good for \n\t"
		<< spGGOp_Analytic_ra->getDataFilename()
		<< " (const iterator accessor, vtemps).";
     else
	 {
	     fail() << spGGOp_Analytic_ra->getDataDescriptor()
		    << " opacity value is out of spec. for \n\t"
		    << spGGOp_Analytic_ra->getDataFilename()
		    << " (const iterator accessor, vtemps).";
	     return false;
	 }


     // ------------------------------------------ //
     // Test the STL-like getOpacity accessor      //
     // const iterator ( density only) for MG data //
     // ------------------------------------------ //

     // clear old data
     for ( int i=0; i<nd*ng; ++i) mgOp[i]=0.0;

     // Calculate the reference solution.
     for ( int id=0; id<nd; ++id )
	 for ( int ig=0; ig<ng; ++ig )
	     vtabulatedOpacity[ id*ng + ig ] 
		 = vdensity[id] * pow( temperature, 4 );
     
     // Obtain new solution
     spGMGOp_Analytic_ra->getOpacity( temperature,
				      cvdensity.begin(),
				      cvdensity.end(), 
				      mgOp.begin() );
     
     // Compare solutions and report the results.
     if ( match( mgOp, vtabulatedOpacity ) )
	 pass() << spGMGOp_Analytic_ra->getDataDescriptor()
		<< " opacity computation was good for \n\t"
		<< spGMGOp_Analytic_ra->getDataFilename()
		<< " (const iterator accessor, vdensity).";
     else
	 {
	     fail() << spGMGOp_Analytic_ra->getDataDescriptor()
		    << " opacity value is out of spec. for \n\t"
		    << spGMGOp_Analytic_ra->getDataFilename()
		    << " (const iterator accessor, vdensity).";
	     return false;
	 }


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

bool tGandolfOpacity::match( const double computedValue,
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

} // end of tGandolfOpacity::match( double, double )

// ------------- //
// Match vectors //
// ------------- //

bool tGandolfOpacity::match( 
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

// ------------------------ //
// Match vectors of vectors //
// ------------------------ //

bool tGandolfOpacity::match( 
    const std::vector< std::vector<double> >& computedValue, 
    const std::vector< std::vector<double> >& referenceValue ) const 
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
} 

} // end namespace rtt_cdi_gandolf_test

//---------------------------------------------------------------------------//
//                            end of tGandolfOpacity.cc
//---------------------------------------------------------------------------//
