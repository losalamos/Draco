//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/GandolfDummyOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:17:16 2001
 * \brief  DummyMultigroupOpacity templated class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DummyMultigroupOpacity.hh"

#include <cmath> // pow(x,n)
#include <vector>

namespace rtt_dummyMultigroupOpacity
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
    DummyMultigroupOpacity::DummyMultigroupOpacity( )
	: dataFilename( "none" ),
  	  dataDescriptor( "DummyMultigroupOpacity" ),
	  energyPolicyDescriptor( "Multigroup" ),
	  numTemperatures( 3 ),
	  numDensities( 2 ),
	  numGroupBoundaries( 4 )
	{
	    temperatureGrid.resize( numTemperatures );
	    densityGrid.resize( numDensities );
	    groupBoundaries.resize( numGroupBoundaries );  
	    
	    for ( int i=0; i<numTemperatures; ++i )
		temperatureGrid[i] = (i+1) * 1.0;
	    for ( int i=0; i<numDensities; ++i )
		densityGrid[i] = (i+1) * 0.1;
	    for ( int i=0; i<numGroupBoundaries; ++i )
		groupBoundaries[i] = 5.0 * pow(10.0,(i-2.0));
	} 
    

    // --------- //
    // Accessors //
    // --------- //

    /*!
     * \brief Opacity accessor that returns a single opacity (or a
     *     vector of opacities for the multigroup EnergyPolicy) that 
     *     corresponds to the provided temperature and density.
     */
std::vector<double> DummyMultigroupOpacity::getOpacity(
    const double targetTemperature,
    const double targetDensity ) const
    { 
	std::vector< double > opacity( numGroupBoundaries-1 );
	for ( int ig=0; ig<numGroupBoundaries-1; ++ig)
	    opacity[ig] = 2.0 * 
		( targetTemperature + targetDensity/1000.0) 
		/ ( groupBoundaries[ig] + groupBoundaries[ig+1] );
	return opacity;
    }

    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     temperatures and a single density value.
     */
std::vector< std::vector<double> > 
DummyMultigroupOpacity::getOpacity(
    const std::vector<double>& targetTemperature,
    const double targetDensity ) const
    { 
	std::vector< std::vector< double > > opacity( targetTemperature.size() );
	
	for ( int it=0; it<targetTemperature.size(); ++it )
	    {
		opacity[it].resize( numGroupBoundaries-1 );
		
		for ( int ig=0; ig<numGroupBoundaries-1; ++ig)
		    opacity[it][ig] = 2.0 * 
			( targetTemperature[it] + targetDensity/1000.0) 
			/ ( groupBoundaries[ig] + groupBoundaries[ig+1] );
	    }

	return opacity;
    }

    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     densities and a single temperature value.
     */
std::vector< std::vector<double> > 
DummyMultigroupOpacity::getOpacity(
    const double targetTemperature,
    const std::vector<double>& targetDensity ) const
    { 
	std::vector< std::vector< double > > opacity( targetDensity.size() );

	for ( int id=0; id<targetDensity.size(); ++id )
	    {
		opacity[id].resize( numGroupBoundaries-1 );
		
		for ( int ig=0; ig<numGroupBoundaries-1; ++ig)
		    opacity[id][ig] = 2.0 * 
			( targetTemperature + targetDensity[id]/1000.0) 
			/ ( groupBoundaries[ig] + groupBoundaries[ig+1] );
	    }

	return opacity;
    }

} // end namespace rtt_dummyMultigroupOpacity


//---------------------------------------------------------------------------//
// end of DummyMultigroupOpacity.cc
//---------------------------------------------------------------------------//
