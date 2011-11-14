//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_ipcress/test/cdi_ipcress_test.hh
 * \author Thomas M. Evans
 * \date   Fri Oct 12 15:36:36 2001
 * \brief  cdi_ipcress test function headers.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_ipcress_test_hh__
#define __cdi_ipcress_test_hh__

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

namespace rtt_cdi_ipcress_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set rtt_cdi_ipcress_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (rtt_cdi_ipcress_test::passed will have its default value of true)

// Needless to say, these can be used in many different combinations or
// ways.  We do not constrain draco tests except that the output must be of
// the form "Test: pass/fail"

bool fail(int line);

bool fail(int line, char *file);

bool pass_msg(const std::string &);

bool fail_msg(const std::string &);

//---------------------------------------------------------------------------//
// PASSING CONDITIONALS
//---------------------------------------------------------------------------//

extern bool passed;

//---------------------------------------------------------------------------//
// DATA EQUIVALENCE FUNCTIONS USED FOR TESTING
//---------------------------------------------------------------------------//

bool match( const double computedValue, const double referenceValue );

//---------------------------------------------------------------------------//

bool match(const std::vector< double >& computedValue, 
	   const std::vector< double >& referenceValue );

//---------------------------------------------------------------------------//

bool match(const std::vector< std::vector<double> >& computedValue, 
	   const std::vector< std::vector<double> >& referenceValue );

//---------------------------------------------------------------------------//

bool match(
	const std::vector< std::vector< std::vector< double > > >& computedValue, 
   const std::vector< std::vector< std::vector< double > > >& referenceValue );
//---------------------------------------------------------------------------//
// COMPARISON FUNCTIONS USED IN IPCRESS OPACITY TESTS
//---------------------------------------------------------------------------//

template < class temperatureType, class densityType, 
	   class testValueType, class opacityClassType >
bool opacityAccessorPassed(const opacityClassType spOpacity, 
			   const temperatureType temperature, 
			   const densityType density, 
			   const testValueType tabulatedValue )
{
    using std::ostringstream;

    // Interpolate the multigroup opacities.
    testValueType grayOpacity
	= spOpacity->getOpacity( temperature, density );
	
    // Make sure that the interpolated value matches previous
    // interpolations. 

    if ( match( grayOpacity, tabulatedValue ) )
    {
	ostringstream message;
	message << spOpacity->getDataDescriptor()
		<< " opacity computation was good for \n\t" 
		<< "\"" << spOpacity->getDataFilename() << "\" data."; 
	pass_msg(message.str());
    }
    else
    {
	ostringstream message;
	message << spOpacity->getDataDescriptor()
		<< " opacity value is out of spec. for \n\t"
		<< "\"" << spOpacity->getDataFilename() << "\" data."; 
	return fail_msg(message.str());
    }
	
    // If we get here then the test passed.
    return true;
}

//---------------------------------------------------------------------------//

template< class opacityClassType >
void testTemperatureGridAccessor(const opacityClassType spOpacity)
{
    using std::ostringstream;

    // Read the temperature grid from the IPCRESS file.     
    std::vector< double > temps = spOpacity->getTemperatureGrid();
	
    // Verify that the size of the temperature grid looks right.  If
    // it is the right size then compare the temperature grid data to 
    // the data specified when we created the IPCRESS file using TOPS.
    if ( temps.size() == spOpacity->getNumTemperatures() &&
	 temps.size() == 3 )
    {
	ostringstream message;
	message << "The number of temperature points found in the data\n\t" 
	        << "grid matches the number returned by the\n\t"
		<< "getNumTemperatures() accessor.";
	pass_msg(message.str());
		
	// The grid specified by TOPS has 3 temperature points.
	std::vector< double > temps_ref( temps.size() );
	temps_ref[0] = 0.1;
	temps_ref[1] = 1.0;
	temps_ref[2] = 10.0;
		
	// Compare the grids.
	if ( match( temps, temps_ref ) )
	{
	    pass_msg("Temperature grid matches.");
	}
	else
	{
	    fail_msg("Temperature grid did not match.");
	}
    }
    else
    {
	ostringstream message;
	message << "The number of temperature points found in the data\n\t"
	        << "grid does not match the number returned by the\n\t"
	        << "getNumTemperatures() accessor. \n"
		<< "Did not test the results returned by\n\t"
		<< "getTemperatureGrid().";
	fail_msg(message.str());
    }
}

//---------------------------------------------------------------------------//
    
template< class opacityClassType >
void testDensityGridAccessor(const opacityClassType spOpacity)
{
    using std::ostringstream;

    // Read the grid from the IPCRESS file.     
    std::vector< double > density = spOpacity->getDensityGrid();
	
    // Verify that the size of the density grid looks right.  If
    // it is the right size then compare the density grid data to 
    // the data specified when we created the IPCRESS file using TOPS.
    if ( density.size() == 3 &&
	 density.size() == spOpacity->getNumDensities() )
    {
	ostringstream message;
	message << "The number of density points found in the data\n\t"
		<< "grid matches the number returned by the\n\t"
		<< "getNumDensities() accessor.";
	pass_msg(message.str());
		
	// The grid specified by TOPS has 3 density points
	std::vector< double > density_ref( density.size() );
	density_ref[0] = 0.1;
	density_ref[1] = 0.5;
	density_ref[2] = 1.0;
		
	// Compare the grids.
	if ( match( density, density_ref ) )
	{
	    pass_msg("Density grid matches.");
	}
	else
	{
	    fail_msg("Density grid did not match.");
	}
    }
    else
    {
	ostringstream message;
	message << "The number of density points found in the data\n\t"
		<< "grid does not match the number returned by the\n\t"
		<< "getNumDensities() accessor. \n"
		<< "Did not test the results returned by\n\t"  
		<< "getDensityGrid().";
	fail_msg(message.str());
    }
}

//---------------------------------------------------------------------------//

template< class opacityClassType >
void testEnergyBoundaryAccessor(const opacityClassType spOpacity)
{
    using std::ostringstream;

    // Read the grid from the IPCRESS file.     
    std::vector< double > ebounds = spOpacity->getGroupBoundaries();

    // Verify that the size of the group boundary grid looks right.  If
    // it is the right size then compare the energy groups grid data to 
    // the data specified when we created the IPCRESS file using TOPS.
    if ( ebounds.size() == 13 &&
	 ebounds.size() == spOpacity->getNumGroupBoundaries() )
    {
	ostringstream message;
	message << "The number of energy boundary points found in the data\n\t"
	        << "grid matches the number returned by the\n\t"
	        << "getNumGroupBoundaries() accessor.";
	pass_msg(message.str());

	// The grid specified by TOPS has 13 energy boundaries.
	std::vector< double > ebounds_ref(ebounds.size());
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
	{
	    pass_msg("Energy group boundary grid matches.");
	}
	else
	{
	    fail_msg("Energy group boundary grid did not match.");
	}    
    }
    else
    {
	ostringstream message;
	message << "The number of energy boundary points found in the data\n\t"
		<< "grid does not match the number returned by the\n\t"
		<< "get NumGroupBoundaries() accessor. \n"
		<< "Did not test the results returned by\n\t"  
		<< "getGroupBoundaries().";
	fail_msg(message.str());
    } 
}

} // end namespace rtt_cdi_ipcress_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS    rtt_cdi_ipcress_test::fail(__LINE__);
#define FAILURE    rtt_cdi_ipcress_test::fail(__LINE__, __FILE__);
#define PASSMSG(a) rtt_cdi_ipcress_test::pass_msg(a);
#define FAILMSG(a) rtt_cdi_ipcress_test::fail_msg(a);

#endif                          // __cdi_ipcress_test_hh__

//---------------------------------------------------------------------------//
//                              end of cdi_ipcress/test/cdi_ipcress_test.hh
//---------------------------------------------------------------------------//
