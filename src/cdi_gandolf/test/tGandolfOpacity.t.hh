//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   tGandolfOpacity.t.hh
 * \author Kelly Thompson
 * \date   Tue Jan 23 14:29:20 2001
 * \brief  Implementation file for tGandolfOpacity
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "tGandolfOpacity.hh"

namespace rtt_cdi_gandolf_test {

// ------------------------- //
// Test the opacity accessor // 
// ------------------------- //

template < class temperatureType, class densityType, 
    class testValueType, class opacityClassType >
bool tGandolfOpacity::opacityAccessorPassed(
    const opacityClassType spOpacity, 
    const temperatureType temperature, 
    const densityType density, 
    const testValueType tabulatedValue )
    {
	// Interpolate the multigroup opacities.
	testValueType grayOpacity
	    = spOpacity->getOpacity( temperature, density );
	
	// Make sure that the interpolated value matches previous
	// interpolations. 

	if ( match( grayOpacity, tabulatedValue ) )
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

// ---------------------------------- //
// Test the Temperature Grid Accessor //
// ---------------------------------- //

template< class opacityClassType >
void tGandolfOpacity::testTemperatureGridAccessor( 
    const opacityClassType spOpacity )
    {
	// Read the temperature grid from the IPCRESS file.     
	std::vector< double > temps = spOpacity->getTemperatureGrid();
	
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
		std::vector< double > temps_ref( temps.size() );
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
    } 

// ------------------------------ //
// Test the Density Grid Accessor //
// ------------------------------ //

template< class opacityClassType >
void tGandolfOpacity::testDensityGridAccessor( 
    const opacityClassType spOpacity )
    {
	// Read the grid from the IPCRESS file.     
	std::vector< double > density = spOpacity->getDensityGrid();
	
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
		std::vector< double > density_ref( density.size() );
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
    } 


// --------------------------------- //
// Test the Energy Boundary Accessor //
// --------------------------------- //

template< class opacityClassType >
void tGandolfOpacity::testEnergyBoundaryAccessor( 
    const opacityClassType spOpacity )
    {
    
     // Read the grid from the IPCRESS file.     
     std::vector< double > ebounds = spOpacity->getGroupBoundaries();

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
// end of tGandolfOpacity.t.hh
//---------------------------------------------------------------------------//
