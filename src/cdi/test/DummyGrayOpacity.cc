//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyGrayOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:33:51 2001
 * \brief  DummyGrayOpacity class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DummyGrayOpacity.hh"
#include <vector>

namespace rtt_dummyGrayOpacity
{

    // ------------ //
    // Constructors //
    // ------------ //

/*!
 * \brief Constructor for DummyOpacity object.
 * 
 * See DummyOpacity.hh for details.
 *
 * Note that everything in this file must be templated by the
 * EnergyPolicy.  All Templated forms of DummyOpacity<EnergyPolicy>
 * must be instantiated in DummyOpacity_pt.cc
 */
    DummyGrayOpacity::DummyGrayOpacity()
	    : dataFilename( "none"), 
	      dataDescriptor( "DummyGrayOpacity" ), 
	      energyPolicyDescriptor( "Gray" ), 
	      numTemperatures( 3 ), 
	      numDensities( 2 )
	{
	    // Set up the temperature and density grid.
	    temperatureGrid.resize( numTemperatures );
	    densityGrid.resize( numDensities );
	    for ( int i=0; i<numTemperatures; ++i )
		temperatureGrid[i] = 1.0 * (i+1);
	    for ( int i=0; i<numDensities; ++i )
		densityGrid[i] = 0.1 * (i+1);
	}

    // --------- //
    // Accessors //
    // --------- //
    
    /*!
     * \brief Opacity accessor that returns a single opacity (or a
     *     vector of opacities for the multigroup EnergyPolicy) that 
     *     corresponds to the provided temperature and density.
     */
double DummyGrayOpacity::getOpacity(
    const double targetTemperature,
    const double targetDensity ) const
    { 
	return targetTemperature + targetDensity/1000.0;
    }
 
std::vector< double > DummyGrayOpacity::getOpacity(
    const std::vector< double >& targetTemperature,
    const double targetDensity ) const
    {
	std::vector< double > grayOpacity(
	    targetTemperature.size() );
	for ( int i=0; i<targetTemperature.size(); ++i )
	    grayOpacity[i] = targetTemperature[i] +
		    targetDensity/1000.0;
	return grayOpacity;
    }
 
std::vector< double > DummyGrayOpacity::getOpacity( 
    const double targetTemperature,
    const std::vector< double >& targetDensity ) const
    {
	std::vector< double > grayOpacity(
	    targetDensity.size() );
	for ( int i=0; i<targetDensity.size(); ++i )
	    grayOpacity[i] = targetTemperature +
		targetDensity[i]/1000.0;
	return grayOpacity;
    }
 
} // end rtt_dummyGrayOpacity

//---------------------------------------------------------------------------//
// end of DummyGrayOpacity.cc
//---------------------------------------------------------------------------//
