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

#include <vector>  // required by match()
#include <cmath>   // define fabs()

#include "tDummyOpacity.hh"
#include "../Release.hh"

#include "DummyGrayOpacity.hh"
#include "DummyGrayOpacity.t.hh"
#include "DummyMultigroupOpacity.hh"
#include "DummyMultigroupOpacity.t.hh"

#include "UnitTestFrame/PassFailStream.hh"
#include "ds++/SP.hh"

// --------------------- //
// Unit Test Frame Stuff //
// --------------------- //

namespace rtt_UnitTestFrame {
    rtt_dsxx::SP<TestApp> TestApp::create( int &argc, char *argv[],
					   std::ostream& os_in ) {
	return rtt_dsxx::SP<TestApp> ( 
	    new rtt_dummy_opacity_test::tDummyOpacity( argc, argv, os_in ));
    }
} // end namespace rtt_UnitTestFrame

// ------------------- //
// tDummyOpacity Stuff //
// ------------------- //

namespace rtt_dummy_opacity_test 
{

tDummyOpacity::tDummyOpacity( int argc, char *argv[], std::ostream& os_in )
    : rtt_UnitTestFrame::TestApp( argc, argv, os_in )
    {
	os() << "Created tDummyOpacity" << std::endl;
    }
 
std::string tDummyOpacity::version() const
{
    return rtt_cdi::release();
}


          //====================================//
          // The DummyOpacity tests start here. //
          //====================================//


std::string tDummyOpacity::runTest()
{

    // ---------------------------- //
    // Create a GrayOpacity object. //
    // ---------------------------- //

    rtt_dsxx::SP< rtt_cdi::GrayOpacity > spDGO;
    
    if ( spDGO = new rtt_dummyGrayOpacity::DummyGrayOpacity() )
	// If we get here then the object was successfully instantiated.
	pass() << "SP to new GrayOpacity object created.";
    else
	{
	    fail() << "Unable to create a SP to new GrayOpacity object.";
	    return "Unable to create a SP to new GrayOpacity object.";
	}

    // ------------------------ //
    // Dummy Gray Opacity Tests //
    // ------------------------ //
    
    double temperature = 0.1; // keV
    double density = 27.0; // g/cm^3
    double tabulatedGrayOpacity = temperature + density/1000.0; // cm^2/g
    
    double opacity = spDGO->getOpacity( temperature, density );

    if ( match ( opacity, tabulatedGrayOpacity ) ) 
	pass() << spDGO->getDataDescriptor()
	       << " getOpacity computation was good.";
    else
	fail() << spDGO->getDataDescriptor()
	       << " getOpacity value is out of spec.";

    // try using a vector of temps.

    std::vector< double > vtemperature(2);
    vtemperature[0] = 0.5; // keV
    vtemperature[1] = 0.7; // keV
    density = 0.35; // g/cm^3
    std::vector< double > vRefOpacity( vtemperature.size() );
    for ( int i=0; i<vtemperature.size(); ++i )
	vRefOpacity[i] = vtemperature[i] + density/1000;

    std::vector< double > vOpacity = spDGO->getOpacity(
	vtemperature, density );

    if ( match ( vOpacity, vRefOpacity ) ) 
	pass() << spDGO->getDataDescriptor()
	       << " getOpacity computation was good for a vector of temps.";
    else
	fail() << spDGO->getDataDescriptor()
	       << " getOpacity value is out of spec. for a vector of temps.";
	

    // try using a vector of densities.

    std::vector< double > vdensity(5);
    vdensity[0] = 0.5;
    vdensity[1] = 1.0;
    vdensity[2] = 3.3;
    vdensity[3] = 5.0;
    vdensity[4] = 27.0;

    vRefOpacity.resize( vdensity.size() );
    for ( int i=0; i<vdensity.size(); ++i )
	vRefOpacity[i] = temperature + vdensity[i]/1000;

    vOpacity = spDGO->getOpacity( temperature, vdensity );

    if ( match ( vOpacity, vRefOpacity ) ) 
	pass() << spDGO->getDataDescriptor()
	       << " getOpacity computation was good for a vector of densities.";
    else
	fail() << spDGO->getDataDescriptor()
	       << " getOpacity value is out of spec. for a vector of densities.";

    // STL-like accessor

    // The virtual base class does not support STL-like accessors, so
    // we must instantiate a DummyOpacity object

    rtt_dsxx::SP< rtt_dummyGrayOpacity::DummyGrayOpacity > spDumGrOp;
    if ( spDumGrOp = new rtt_dummyGrayOpacity::DummyGrayOpacity() )
	pass() << "SP to new DummyGrayOpacity object created.";
    else
	{
	    fail() << "Unable to create a SP to a new DummyGrayOpacity object.";
	    return "Unable to create a SP to a new DummyGrayOpacity object.";
	}

    vOpacity.resize( vdensity.size() );

    spDumGrOp->getOpacity( temperature, vdensity.begin(), vdensity.end(),
			   vOpacity.begin() );

    

    // ----------------------------------------- //
    // Create a Dummy Multigroup Opacity object. //
    // ----------------------------------------- //

    rtt_dsxx::SP< rtt_cdi::MultigroupOpacity > spDmgO;
    
    if ( spDmgO = new rtt_dummyMultigroupOpacity::DummyMultigroupOpacity() )
	// If we get here then the object was successfully instantiated.
	pass() << "SP to new MultigroupOpacity object created.";

    // --------------- //
    // MG Opacity test //
    // --------------- //
    
    // Setup the test point.
    temperature = 0.01; // keV
    density = 2.0; // g/cm^3

    // The dummy opacity object should have 3 groups.  Check it.
    int ng = spDmgO->getNumGroupBoundaries()-1;
    if ( ng == 3 )
	pass() << "Correct number of groups found for MultigroupOpacity object.";
    else
	fail() << "Wrong number of groups found for MultigroupOpacity object.";

    const std::vector< double > energyBoundaries =
	spDmgO->getGroupBoundaries();

    // Create a container that hold all the MG opacities for a
    // specified temperature and density.  Fill this container with
    // the values that DummyMultigroupOpacity should contain.
    std::vector< double > tabulatedMGOpacity( ng );
    for ( int ig=0; ig<ng; ++ig)
	tabulatedMGOpacity[ig] = 2*(temperature + density/1000)
	    / (energyBoundaries[ig]+energyBoundaries[ig+1]);

    // Use the getOpacity accessor to obtain the MG opacities for a
    // specified temperature and density.
    std::vector< double > mgOpacity =
	spDmgO->getOpacity( temperature, density );

    // Make sure the accessor values match the expected values.
    if ( match( mgOpacity, tabulatedMGOpacity ) )
	pass() << spDmgO->getDataDescriptor()
	       << " getOpacity computation was good.";
    else
	fail() << spDmgO->getDataDescriptor()
	       << " getOpacity value is out of spec.";
    
    // Repeat with a vector of temps.
    
    // Reference values.

    // The opacity container is a vector<vector<double>>.  Each nested 
    // vector contains all of the group opacity values for a single
    // temperature. 

    // a MG opacity set for a single temperature, density combination
    // can be extracted from this container by using the following
    // type of assignment.
    // std::vector< double > vec1 = vRefMgOpacity[0];

    std::vector< std::vector< double > > vRefMgOpacity( ng );
    for ( int it=0; it<vtemperature.size(); ++it )
	{
	    vRefMgOpacity[it].resize( ng );
	    for ( int ig=0; ig<ng; ++ig )
		vRefMgOpacity[it][ig] = 
		    2.0 * ( vtemperature[it] + density/1000.0 ) 
		    / ( energyBoundaries[ig] + energyBoundaries[ig+1] );
	}
    
    // Retrieve the same set of opacity values via the getOpacity() accessor.
    std::vector< std::vector< double > > vMgOpacity 
	= spDmgO->getOpacity( vtemperature, density );
    
    // Compare the results.
    if ( match( vMgOpacity, vRefMgOpacity ) )
	pass() << spDmgO->getDataDescriptor()
	       << " getOpacity computation was good for a vector of  temps.";
    else
	fail() << spDmgO->getDataDescriptor()
	       << " getOpacity value is out of spec. for a vector of temps.";
    

    // STL-like accessor (MG opacities)

    // The STL accessors are only available to dummyOpacity objects
    // (not generic opacity objects).

    rtt_dsxx::SP< rtt_dummyMultigroupOpacity::DummyMultigroupOpacity > spDumMgOp;
    if ( spDumMgOp = new rtt_dummyMultigroupOpacity::DummyMultigroupOpacity() )
	pass() << "SP to new DummyMultigroupOpacity object created.";
    else
	{
	    fail() << "Unable to create a SP to a new DummyMultigroupOpacity object.";
	    return "Unable to create a SP to a new DummyMultigroupOpacity object.";
	}

    // The STL-like accessors only work with 1-D containers.
    vOpacity.resize( vtemperature.size() * ng );
    vRefOpacity.resize( vtemperature.size() * ng );

    // vdensity.size() == vtemperature.size()
    vdensity.resize( vtemperature.size() );
    
    // Reference Values
    for ( int it=0; it<vtemperature.size(); ++it )
	for ( int ig=0; ig<ng; ++ig )
	    vRefOpacity[it*ng+ig] = 
		2.0 * ( vtemperature[it] + vdensity[it]/1000.0 )
		/ ( energyBoundaries[ig] + energyBoundaries[ig+1] );

    // Obtain values using getOpacity() accessor.
    spDumMgOp->getOpacity( vtemperature.begin(), vtemperature.end(),
			   vdensity.begin(), vdensity.end(),
			   vOpacity.begin() );
    
    // Compare the results:
    if ( match( vOpacity, vRefOpacity ) )
    	pass() << spDumMgOp->getDataDescriptor()
	       << " STL getOpacity() computation was good for a\n"
	       << " vector of temps. and a vector of densities.";
    else
	fail() << spDumMgOp->getDataDescriptor()
	       << " STL getOpacity() value is out of spec. for a\n"
	       << " vector of temps. and a vector of densities.";


//     std::cout << "MG STL accessor:" << std::endl << std::endl;
//     for ( int it=0; it<vtemperature.size(); ++it )
// 	{
// 	    std::cout << "Temperature = " << vtemperature[it] << " keV." << std::endl
// 		      << "Density     = " << vdensity[it]     << " g/cm^3." << std::endl;
// 	    for ( int ig=0; ig<ng; ++ig )
// 		{
// 		    refOpacity = 2.0 * ( vtemperature[it] + vdensity[it]/1000.0 )
// 			/ ( energyBoundaries[ig] + energyBoundaries[ig+1] );
// 		    std::cout << "   ig = " << ig 
// 			      << "   Opacity = " 
// 			      << vOpacity[it*ng+ig] << " cm^2/g."
// 			      << "   RefOpacity = "
// 			      << refOpacity << std::endl;
// 		}
// 	}


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

bool tDummyOpacity::match( double computedValue,
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

} // end of tDummyOpacity::match( double, double )

bool tDummyOpacity::match( const std::vector< double >& computedValue, 
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
} // end of tDummyOpacity::match( vector<double>, vector<double> )

bool tDummyOpacity::match( 
    const std::vector< std::vector< double > >& computedValue, 
    const std::vector< std::vector< double > >& referenceValue ) const 
{
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
		    // If the comparison fails then stop testing and
		    // return "false" to indicate that the test
		    // failed. 
		    if ( reldiff > TOL ) return false; 
		}
	}
    return true;
} // end of tDummyOpacity::match( vector<double>, vector<double> )

} // end namespace rtt_dummy_opacity_test

//---------------------------------------------------------------------------//
// end of tDummyOpacity.cc
//---------------------------------------------------------------------------//
