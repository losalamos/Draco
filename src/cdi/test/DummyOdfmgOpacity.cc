//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyOdfmgOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:17:16 2001
 * \brief  DummyOdfmgOpacity templated class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DummyOdfmgOpacity.hh"

#include <cmath> // pow(x,n)
#include "ds++/Assert.hh"

namespace rtt_cdi_test
{

// ------------ //
// Constructors //
// ------------ //

/*!
 * \brief Constructor for DummyOdfmgOpacity object.
 * 
 * \sa The constructor assigns fixed values for all of the member
 *     data.  Every instance of this object has the same member
 *     data. 
 *
 *     Temperatures     = { 1.0, 2.0, 3.0 }
 *     Densities        = { 0.1, 0.2 }
 *     EnergyBoundaries = { 0.05, 0.5, 5.0, 50.0 }
 *     BandBoundaries   = { 0.00, 0.125, 0.25, 0.50, 1.00 }
 *
 */
DummyOdfmgOpacity::DummyOdfmgOpacity(rtt_cdi::Reaction reaction,
		rtt_cdi::Model    model)
	: dataFilename( "none" ),
	dataDescriptor( "DummyOdfmgOpacity" ),
	energyPolicyDescriptor( "Odfmg" ),
	numTemperatures( 3 ),
	numDensities( 2 ),
	numGroupBoundaries( 4 ),
	numBandBoundaries( 5 ),
	reaction_type(reaction),
	model_type(model)
{
	temperatureGrid.resize( numTemperatures );
	densityGrid.resize( numDensities );
	groupBoundaries.resize( numGroupBoundaries );
	bandBoundaries.resize( numBandBoundaries );

	for ( int i=0; i<numTemperatures; ++i )
		temperatureGrid[i] = (i+1) * 1.0;
	for ( int i=0; i<numDensities; ++i )
		densityGrid[i] = (i+1) * 0.1;

	// due to the way the CDI test is implemented, these group boundaries MUST
	// match the multigroup frequency boundaries
	for ( int i=0; i<numGroupBoundaries; ++i )
		groupBoundaries[i] = 5.0 * std::pow(10.0,(i-2.0));

	bandBoundaries[0] = 0.0;
	for ( int i=1; i < numBandBoundaries; ++i )
		bandBoundaries[i] = std::pow(2,(i - 4.0));
	//Check(bandBoundaries[numBandBoundaries - 1] == 1.0);
} 

// Constructor for entering a different group boundary structure than the
// default 
DummyOdfmgOpacity::DummyOdfmgOpacity(
		rtt_cdi::Reaction reaction,
		rtt_cdi::Model    model,
		int num_groupboundaries, int num_bandboundaries)
	: dataFilename( "none" ),
	dataDescriptor( "DummyOdfmgOpacity" ),
	energyPolicyDescriptor( "Odfmg" ),
	numTemperatures( 3 ),
	numDensities( 2 ),
	numGroupBoundaries( num_groupboundaries ),
	numBandBoundaries( num_bandboundaries ),
	reaction_type(reaction),
	model_type(model)
{
	temperatureGrid.resize( numTemperatures );
	densityGrid.resize( numDensities );
	groupBoundaries.resize( numGroupBoundaries );  
	bandBoundaries.resize( numBandBoundaries );

	for ( int i=0; i<numTemperatures; ++i )
		temperatureGrid[i] = (i+1) * 1.0;
	for ( int i=0; i<numDensities; ++i )
		densityGrid[i] = (i+1) * 0.1;
	for ( int i=0; i<numGroupBoundaries; ++i )
		groupBoundaries[i] = 5.0 * std::pow(10.0,(i-2.0));

	bandBoundaries[0] = 0.0;
	for ( int i=1; i < numBandBoundaries; ++i )
		bandBoundaries[i] = std::pow(2,(i - numBandBoundaries - 1.0));
	//Check(bandBoundaries[numBandBoundaries - 1] == 1.0);
} 

// --------- //
// Accessors //
// --------- //

/*!
 * \brief Opacity accessor that returns a vector of opacities (one 
 *     for each group) that corresponds to the provided
 *     temperature and density.   
 *
 *     Opacity = 2 * ( temperature + density/1000 ) 
 *                 / ( E_high + E_low ) * 10^(band - 2)
 *
 */
std::vector< std::vector<double> > DummyOdfmgOpacity::getOpacity( 
		double targetTemperature,
		double targetDensity ) const
{ 
	const int numGroups = getNumGroups();
	const int numBands  = getNumBands();

	std::vector< std::vector<double> > opacity( numGroups );
	for ( int group = 0; group < numGroups; group++)
	{
		opacity[group].resize(numBands);

		for (int band = 0; band < numBands; band++)
		{
			opacity[group][band] = 2.0 * 
				( targetTemperature + targetDensity/1000.0) 
				/ ( groupBoundaries[group] + groupBoundaries[group+1] )
				* pow(10.0, band - 2);
		}
	}
	return opacity;
}

/*!
 * \brief Opacity accessor that returns a vector of multigroup
 *     opacities corresponding to the provided vector of
 *     temperatures and a single density.  Each multigroup opacity 
 *     is in itself a vector of numGroups opacities.
 *
 *     Opacity = 2 * ( temperature + density/1000 ) 
 *                 / ( E_high + E_low )
 */
std::vector< std::vector< std::vector<double> > > 
DummyOdfmgOpacity::getOpacity( 
		const std::vector<double>& targetTemperature,
		double targetDensity ) const
{ 
	std::vector< std::vector< std::vector<double> > > opacity( targetTemperature.size() );

	for ( int i=0; i<targetTemperature.size(); ++i )
	{
		opacity[i] = getOpacity(targetTemperature[i], targetDensity);
	}
	return opacity;
}

/*!
 * \brief Opacity accessor that returns a vector of multigroup
 *     opacities corresponding to the provided vector of
 *     densities and a single temperature.  Each multigroup opacity 
 *     is in itself a vector of numGroups opacities.
 *
 *     Opacity = 2 * ( temperature + density/1000 ) 
 *                 / ( E_high + E_low )
 */
std::vector< std::vector< std::vector<double> > > 
DummyOdfmgOpacity::getOpacity(
		double targetTemperature,
		const std::vector<double>& targetDensity ) const
{ 
	std::vector< std::vector< std::vector<double> > > opacity( targetDensity.size() );

	//call our regular getOpacity function for every target density
	for ( int i=0; i<targetDensity.size(); ++i )
	{
		opacity[i] = getOpacity(targetTemperature, targetDensity[i]);
	}
	return opacity;
}

} // end namespace rtt_cdi_test

//---------------------------------------------------------------------------//
// end of DummyOdfmgOpacity.cc
//---------------------------------------------------------------------------//
