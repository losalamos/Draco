//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfOpacity.cc
 * \author Kelly Thompson
 * \date   Wed Jul 12 16:11:55 2000
 * \brief  GandolfOpacity templated class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfOpacity.hh"

#include "GandolfWrapper.hh"    // we make calls to the wrapper routines.

#include "GandolfFile.hh"       // we have smart pointers to
                                // GandolfFile objects.

#include "GandolfException.hh"  // Since we call wrapper routines we
                                // need to be able to throw exceptions
                                // if the Gandolf libraries return an 
                                // error.

#include "GandolfDataTable.hh"  // we have a smart pointer to a
                                // GandolfDataTable object.

#include <cmath> // we need to define log(double) and exp(double)

#include "ds++/Assert.hh" // we make use of Require()

namespace rtt_cdi_gandolf
{

    // ------------ //
    // Constructors //
    // ------------ //

/*!
 * \brief Constructor for GandolfOpacity object.
 * 
 * See GandolfOpacity.hh for details.
 *
 * Note that everything in this file must be templated by the
 * EnergyPolicy.  All Templated forms of GandolfOpacity<EnergyPolicy>
 * must be instantiated in GandolfOpacity_pt.cc
 */
template < class EnergyPolicy >
GandolfOpacity< EnergyPolicy >::GandolfOpacity( 
    const rtt_dsxx::SP<GandolfFile> _spGandolfFile,
    const int _materialID,
    const Model _opacityModel,
    const Reaction _opacityReaction )
    : spGandolfFile( _spGandolfFile ),
    materialID( _materialID ),
    numKeys( 0 ),
    opacityModel( _opacityModel ),
    opacityReaction( _opacityReaction )
    {
	// Verify that the requested material ID is available in the
	// specified IPCRESS file.
	if ( ! spGandolfFile->materialFound( materialID ) )
	    throw gkeysException( -1 );

	// Retrieve keys available fo this material from the IPCRESS
	// file.  wgkeys() returns vKnownKeys, numKeys and errorCode.
	int errorCode;
	wrapper::wgkeys( spGandolfFile->getDataFilename(), materialID, 
			 vKnownKeys, wrapper::maxKeys, numKeys,
			 errorCode );
	if ( errorCode !=0 ) throw gkeysException( errorCode );

	// Create the data table object and fill it with the table
	// data from the IPCRESS file.
 	spGandolfDataTable = new GandolfDataTable(
 	    getEnergyPolicyDescriptor(),
	    opacityModel, 
	    opacityReaction,
	    vKnownKeys,
	    materialID, 
	    spGandolfFile );

    } // end of GandolfData constructor
 
    /*!
     * \brief Default GandolfOpacity() destructor.
     *
     * This is required to correctly release memory when a
     * GandolfOpacity<EnergyPolicy> is destroyed.
     */
template < class EnergyPolicy >
GandolfOpacity< EnergyPolicy >::~GandolfOpacity()
    {
	 // empty
    }


    // --------- //
    // Accessors //
    // --------- //

    /*!
     * \brief Returns a "plain English" description of the opacity
     *     data that this class references. (e.g. "Gray Rosseland
     *     Scattering".) 
     *
     * The definition of this function is not included here to prevent 
     *     the inclusion of the GandolfFile.hh definitions within this 
     *     header file.
     */
template < class EnergyPolicy >
const std::string& GandolfOpacity< EnergyPolicy >::getDataDescriptor() const 
    {
	// call the correct function from the GandolfDataTable
	// object. 
	return spGandolfDataTable->getDataDescriptor(); 
    }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of (temperature,density) tuples.
     *     An opacity value will be returned for each tuple.  The
     *     temperatureIterator and density iterators are required to
     *     be the same length.  The opacity iterator should also have
     *     this same length for gray data or this length times the
     *     number of energy groups for multigroup data.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     */
template < class EnergyPolicy >
template < class InputIterator, class OutputIterator >
OutputIterator 
GandolfOpacity< EnergyPolicy >::getOpacity(
    InputIterator temperatureIterator, 
    InputIterator temperatureIteratorEnd,
    InputIterator densityIterator, 
    InputIterator densityIteratorEnd,
    OutputIterator opacityIterator ) const
    { 
	// call the correct accessor from the PolicyTemplate
	return EnergyPolicy::getOpacity(
	    temperatureIterator, temperatureIteratorEnd,
	    densityIterator, densityIteratorEnd,
	    opacityIterator, spGandolfDataTable );
    }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of temperatures in an STL container.
     *     An opacity value will be returned for each temperature
     *     provided.  The opacity iterator should be the same length
     *     as the temperatureIterator for gray data or the length of
     *     the temperatureIterator times the number of energy groups
     *     for multigroup data.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     */
template < class EnergyPolicy >
template < class InputIterator, class OutputIterator >
OutputIterator 
GandolfOpacity< EnergyPolicy >::getOpacity(
    InputIterator temperatureIterator, 
    InputIterator temperatureIteratorEnd,
    const double targetDensity,
    OutputIterator opacityIterator ) const
    { 
	// call the correct accessor from the PolicyTemplate
	return EnergyPolicy::getOpacity(
	    temperatureIterator, temperatureIteratorEnd,
	    targetDensity,
	    opacityIterator, spGandolfDataTable );
    }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of densities in an STL container.
     *     An opacity value will be returned for each density
     *     provided.  The opacity iterator should be the same length
     *     as the densityIterator for gray data or the length of the
     *     densityIterator times the number of energy groups for
     *     multigroup data.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     */
template < class EnergyPolicy >
template < class InputIterator, class OutputIterator >
OutputIterator 
GandolfOpacity< EnergyPolicy >::getOpacity(
    const double targetTemperature,
    InputIterator densityIterator, 
    InputIterator densityIteratorEnd,
    OutputIterator opacityIterator ) const
    { 
	// call the correct accessor from the PolicyTemplate
	return EnergyPolicy::getOpacity(
	    targetTemperature,
	    densityIterator, densityIteratorEnd,
	    opacityIterator, spGandolfDataTable );
    }

    /*!
     * \brief Opacity accessor that returns a single opacity (or a
     *     vector of opacities for the multigroup EnergyPolicy) that 
     *     corresponds to the provided temperature and density.
     */
template < class EnergyPolicy >
GandolfOpacity< EnergyPolicy >::OpacityType 
GandolfOpacity< EnergyPolicy >::getOpacity(
    const double targetTemperature,
    const double targetDensity ) const
    { 
	// call the correct accessor from the PolicyTemplate
	return EnergyPolicy::getOpacity( targetTemperature,
					 targetDensity,
					 spGandolfDataTable ); 
    }

    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     temperatures and a single density value.
     */
template < class EnergyPolicy >
std::vector< typename GandolfOpacity< EnergyPolicy >::OpacityType > 
GandolfOpacity< EnergyPolicy >::getOpacity(
    const std::vector<double>& targetTemperature,
    const double targetDensity ) const
    { 
	// call the correct accessor from the PolicyTemplate
	return EnergyPolicy::getOpacity( targetTemperature,
					 targetDensity,
					 spGandolfDataTable ); 
    }

    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     densities and a single temperature value.
     */
template < class EnergyPolicy >
std::vector< typename GandolfOpacity< EnergyPolicy >::OpacityType > 
GandolfOpacity< EnergyPolicy >::getOpacity(
    const double targetTemperature,
    const std::vector<double>& targetDensity ) const
    { 
	// call the correct accessor from the PolicyTemplate
	return EnergyPolicy::getOpacity( targetTemperature,
					 targetDensity,
					 spGandolfDataTable ); 
    }

    // It is not clear how to assume order of opacity(temp,dens) when
    // accessed in this manner --> for now use the STL-style accessor
    // or a loop over one of the other vector-accessors.

// template < class EnergyPolicy >
// std::vector< typename GandolfOpacity< EnergyPolicy >::OpacityType > 
// GandolfOpacity< EnergyPolicy >::getOpacity(
//     const std::vector<double>& targetTemperature,
//     const std::vector<double>& targetDensity ) const
//     { 
// 	return EnergyPolicy::getOpacity( targetTemperature,
// 					 targetDensity,
// 					 spGandolfDataTable ); 
//     }

    /*!
     * \brief Returns a vector of temperatures that define the cached
     *     opacity data table.
     */
template < class EnergyPolicy >
std::vector<double> GandolfOpacity< EnergyPolicy >::getTemperatureGrid() const
    {
	// Retrieve the temperature grid for the Data Object.  These
	// are log(temp) values.
	std::vector<double> tGrid =
	    spGandolfDataTable->getLogTemperatures();
	// convert the log(temp) grid to a temp grid.
	for ( int i=0; i<tGrid.size(); ++i )
	    tGrid[i] = exp( tGrid[i] );
	// return the temperature grid.
	return tGrid;
    }

    /*!
     * \brief Returns the size of the temperature grid.
     */
template < class EnergyPolicy >
int GandolfOpacity< EnergyPolicy >::getNumTemperatures() const
    {
	return spGandolfDataTable->getNumTemperatures();
    }

    /*!
     * \brief Returns a vector of densities that define the cached
     *     opacity data table.
     */
template < class EnergyPolicy >
std::vector<double> GandolfOpacity< EnergyPolicy >::getDensityGrid() const
    {
	// Retrieve the density grid for the Data Object.  These
	// are log(rho) values.
	std::vector<double> rhoGrid =
	    spGandolfDataTable->getLogDensities();
	// convert the log(rho) grid to a density grid.
	for ( int i=0; i<rhoGrid.size(); ++i )
	    rhoGrid[i] = exp( rhoGrid[i] );
	// return the density grid.
	return rhoGrid;
    }

    /*! 
     * \brief Returns the size of the density grid.
     */
template < class EnergyPolicy >
int GandolfOpacity< EnergyPolicy >::getNumDensities() const
    {
	return spGandolfDataTable->getNumDensities();
    }

    /*!
     * \brief Returns a vector of energy values (keV) that define the
     *     energy boundaries of the cached multigroup opacity data
     *     table.  (This accessor is only valid for the Multigroup
     *     EnergyPolicy version of GandolfOpacity.)
     */
template < class EnergyPolicy >
const std::vector<double>& GandolfOpacity< EnergyPolicy >::getGroupBoundaries() const
    {
	return spGandolfDataTable->getGroupBoundaries();
    }

    /*!
     * \brief Returns the number of group boundaries found in the
     *     current multigroup data set.
     */
template < class EnergyPolicy >
int GandolfOpacity< EnergyPolicy >::getNumGroupBoundaries() const
    {
	return spGandolfDataTable->getNumGroupBoundaries();
    }

		       // ----------------- //
		       // Gray Policy Class //
		       // ----------------- //

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of (temperature,density) tuples.
     *     A set of gray opacity values will be returned for each
     *     tuple.  The temperatureIterator and densityIterator
     *     are required to be the same length.  The opacityIterator
     *     should have a length equal to the the temperatureIterator.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     */
template < class InputIterator, class OutputIterator >
OutputIterator Gray::getOpacity(
    InputIterator temperatureIterator, 
    InputIterator temperatureIteratorEnd,
    InputIterator densityIterator, 
    InputIterator densityIteratorEnd,
    OutputIterator opacityIterator,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	// from twix:/scratch/tme/kai/KCC_BASE/include/algorithm

	// assert that the two input iterators have compatible sizes.
 	Require( std::distance( temperatureIterator, temperatureIteratorEnd )
 		 == std::distance( densityIterator, densityIteratorEnd ) );

	// Loop over all (temperature,density) tuple values.
	for ( ; temperatureIterator != temperatureIteratorEnd;
	      ++temperatureIterator, ++densityIterator, ++opacityIterator )
	    // Call the Gandolf Logorithmic Interpolator.
	    wrapper::wgintgrlog( spGandolfDataTable->getLogTemperatures(),
				 spGandolfDataTable->getNumTemperatures(), 
				 spGandolfDataTable->getLogDensities(), 
				 spGandolfDataTable->getNumDensities(),
				 spGandolfDataTable->getLogOpacities(), 
				 spGandolfDataTable->getNumOpacities(),
				 log(*temperatureIterator),
				 log(*densityIterator), 
				 *opacityIterator );
	return opacityIterator;
    }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of temperatures in an STL container.
     *     An set of opacity values will be returned for each temperature
     *     provided.  The opacity iterator should have length equal to 
     *     the length of the temperatureIterator.
     */
template < class InputIterator, class OutputIterator >
OutputIterator Gray::getOpacity(
    InputIterator temperatureIterator, 
    InputIterator temperatureIteratorEnd,
    const double targetDensity,
    OutputIterator opacityIterator,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	for ( ; temperatureIterator != temperatureIteratorEnd;
	      ++temperatureIterator, ++opacityIterator )
	    wrapper::wgintgrlog( spGandolfDataTable->getLogTemperatures(),
				 spGandolfDataTable->getNumTemperatures(), 
				 spGandolfDataTable->getLogDensities(), 
				 spGandolfDataTable->getNumDensities(),
				 spGandolfDataTable->getLogOpacities(), 
				 spGandolfDataTable->getNumOpacities(),
				 log( *temperatureIterator ),
				 log( targetDensity ), 
				 *opacityIterator );
	return opacityIterator;
    }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of densities in an STL container.
     *     An opacity value will be returned for each density
     *     provided.  The opacity iterator should have length equal to 
     *     to the length of the densityIterator.
     */
template < class InputIterator, class OutputIterator >
OutputIterator Gray::getOpacity(
    const double targetTemperature,
    InputIterator densityIterator, 
    InputIterator densityIteratorEnd,
    OutputIterator opacityIterator,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	for ( ; densityIterator != densityIteratorEnd;
	      ++densityIterator, ++opacityIterator )
	    wrapper::wgintgrlog( spGandolfDataTable->getLogTemperatures(),
				 spGandolfDataTable->getNumTemperatures(), 
				 spGandolfDataTable->getLogDensities(), 
				 spGandolfDataTable->getNumDensities(),
				 spGandolfDataTable->getLogOpacities(), 
				 spGandolfDataTable->getNumOpacities(),
				 log( targetTemperature ),
				 log( *densityIterator ), 
				 *opacityIterator );
	return opacityIterator;
    }

   /*!
     * \brief Opacity accessor that returns a single opacity that 
     *     corresponds to the provided temperature and density.
     */
Gray::OpacityType Gray::getOpacity( 
    const double targetTemperature,
    const double targetDensity,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	OpacityType opacity;
	// logorithmic interpolation:
	wrapper::wgintgrlog( spGandolfDataTable->getLogTemperatures(),
			     spGandolfDataTable->getNumTemperatures(), 
			     spGandolfDataTable->getLogDensities(), 
			     spGandolfDataTable->getNumDensities(),
			     spGandolfDataTable->getLogOpacities(), 
			     spGandolfDataTable->getNumOpacities(),
			     log(targetTemperature),
			     log(targetDensity), 
			     opacity );
	return opacity;
    }

    /*!
     * \brief Opacity accessor that returns a vector of
     *     opacities that correspond to the provided vector of
     *     temperatures and a single density value.
     */
std::vector< Gray::OpacityType > Gray::getOpacity( 
    const std::vector<double>& targetTemperature,
    const double targetDensity,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	std::vector< OpacityType > opacity( targetTemperature.size() );
	for ( int i=0; i<targetTemperature.size(); ++i )
	    // logorithmic interpolation:
	    wrapper::wgintgrlog( spGandolfDataTable->getLogTemperatures(),
				 spGandolfDataTable->getNumTemperatures(), 
				 spGandolfDataTable->getLogDensities(), 
				 spGandolfDataTable->getNumDensities(),
				 spGandolfDataTable->getLogOpacities(), 
				 spGandolfDataTable->getNumOpacities(),
				 log(targetTemperature[i]),
				 log(targetDensity), 
				 opacity[i] );
	return opacity;
    }

    /*!
     * \brief Opacity accessor that returns a vector of
     *     opacities that correspond to the provided vector of
     *     densities and a single temperature value.
     */
std::vector< Gray::OpacityType > Gray::getOpacity( 
    const double targetTemperature,
    const std::vector<double>& targetDensity,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	std::vector< OpacityType > opacity( targetDensity.size() );
	for ( int i=0; i<targetDensity.size(); ++i )
	    // logorithmic interpolation:
	    wrapper::wgintgrlog( spGandolfDataTable->getLogTemperatures(),
				 spGandolfDataTable->getNumTemperatures(), 
				 spGandolfDataTable->getLogDensities(), 
				 spGandolfDataTable->getNumDensities(),
				 spGandolfDataTable->getLogOpacities(), 
				 spGandolfDataTable->getNumOpacities(),
				 log(targetTemperature),
				 log(targetDensity[i]), 
				 opacity[i] );
	return opacity;
    }

    // This routine was removed because of confusion surrounding the
    // ordering of temperature, density and multigroup indicies in the 
    // returned opacity container.

// std::vector< Gray::OpacityType > Gray::getOpacity( 
//     const std::vector<double>& targetTemperature,
//     const std::vector<double>& targetDensity,
//     const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
//     {
// 	int nt = targetTemperature.size();
// 	int nd = targetDensity.size();
// 	std::vector< OpacityType > opacity( nt*nd );
// 	for ( int i=0; i<nt; ++i )
// 	    for ( int j=0; j<nd; ++j )
// 		// logorithmic interpolation:
// 		wrapper::wgintgrlog( spGandolfDataTable->getLogTemperatures(),
// 				     spGandolfDataTable->getNumTemperatures(), 
// 				     spGandolfDataTable->getLogDensities(), 
// 				     spGandolfDataTable->getNumDensities(),
// 				     spGandolfDataTable->getLogOpacities(), 
// 				     spGandolfDataTable->getNumOpacities(),
// 				     log(targetTemperature[i]),
// 				     log(targetDensity[j]), 
// 				     opacity[i*nd+j] );
// 	return opacity;
//     }



		    // ----------------------- //
		    // Multigroup Policy Class //
		    // ----------------------- //

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of (temperature,density) tuples.
     *     A set of multigroup opacity values will be returned for
     *     each tuple.  The temperatureIterator and density iterators
     *     are required to be the same length.  The opacity iterator
     *     should have a length equal to the the temperatureIterator
     *     multiplied by the number of energy groups.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     */
template < class InputIterator, class OutputIterator >
OutputIterator Multigroup::getOpacity(
    InputIterator temperatureIterator, 
    InputIterator temperatureIteratorEnd,
    InputIterator densityIterator, 
    InputIterator densityIteratorEnd,
    OutputIterator opacityIterator,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	// from twix:/scratch/tme/kai/KCC_BASE/include/algorithm

	// assert that the two input iterators have compatible sizes.
 	Require( std::distance( temperatureIterator, temperatureIteratorEnd )
 		 == std::distance( densityIterator, densityIteratorEnd ) );

	// number of groups in this multigroup set.
	const int ng = spGandolfDataTable->getNumGroupBoundaries()-1;
	
	// temporary opacity vector used by the wrapper.  The returned 
	// data will be copied into the opacityIterator.
	std::vector<double> opacity(ng);

	// loop over the (temperature,density) tuple.
	for ( ; temperatureIterator != temperatureIteratorEnd;
	      ++temperatureIterator, ++densityIterator )
	    {
		// Call the Gandolf interpolator.
		// the vector opacity is returned.
		wrapper::wgintmglog( spGandolfDataTable->getLogTemperatures(),
				     spGandolfDataTable->getNumTemperatures(), 
				     spGandolfDataTable->getLogDensities(), 
				     spGandolfDataTable->getNumDensities(),
				     spGandolfDataTable->getNumGroupBoundaries(),
				     spGandolfDataTable->getLogOpacities(), 
				     spGandolfDataTable->getNumOpacities(),
				     log(*temperatureIterator),
				     log(*densityIterator), 
				     opacity );
		
		// The opacity vector contains the solution.  Now
		// we copy this solution into the OutputIterator
		for ( int i=0; i<ng; ++i, ++opacityIterator )
		    *opacityIterator = opacity[i];
	    }
	return opacityIterator;
    }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of temperatures in an STL container.
     *     An set of opacity values will be returned for each temperature
     *     provided.  The opacity iterator should have length equal to 
     *     the length of the temperatureIterator times the number of
     *     energy groups.
     */
template < class InputIterator, class OutputIterator >
OutputIterator Multigroup::getOpacity(
    const double targetTemperature,
    InputIterator densityIterator, 
    InputIterator densityIteratorEnd,
    OutputIterator opacityIterator,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	for ( ; densityIterator != densityIteratorEnd; ++densityIterator )
	    {
		const int ng = spGandolfDataTable->getNumGroupBoundaries()-1;
		std::vector<double> opacity(ng);
		wrapper::wgintmglog( spGandolfDataTable->getLogTemperatures(),
				     spGandolfDataTable->getNumTemperatures(), 
				     spGandolfDataTable->getLogDensities(), 
				     spGandolfDataTable->getNumDensities(),
				     spGandolfDataTable->getNumGroupBoundaries(),
				     spGandolfDataTable->getLogOpacities(), 
				     spGandolfDataTable->getNumOpacities(),
				     log( targetTemperature ),
				     log( *densityIterator ), 
				     opacity );
		
		// The opacity vector contains the solution.  Now
		// we copy this solution into the OutputIterator
		for ( int i=0; i<ng; ++i, ++opacityIterator )
		    *opacityIterator = opacity[i];
		
	    }
	return opacityIterator;
    }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of densities in an STL container.
     *     An opacity value will be returned for each density
     *     provided.  The opacity iterator should have length equal to 
     *     to the length of the densityIterator times the number of
     *     energy groups.
     */
template < class InputIterator, class OutputIterator >
OutputIterator Multigroup::getOpacity(
    InputIterator temperatureIterator, 
    InputIterator temperatureIteratorEnd,
    const double targetDensity,
    OutputIterator opacityIterator,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	for ( ; temperatureIterator != temperatureIteratorEnd;
	     ++temperatureIterator )
	    {
		const int ng = spGandolfDataTable->getNumGroupBoundaries()-1;
		std::vector<double> opacity(ng);
		wrapper::wgintmglog( spGandolfDataTable->getLogTemperatures(),
				     spGandolfDataTable->getNumTemperatures(), 
				     spGandolfDataTable->getLogDensities(), 
				     spGandolfDataTable->getNumDensities(),
				     spGandolfDataTable->getNumGroupBoundaries(),
				     spGandolfDataTable->getLogOpacities(), 
				     spGandolfDataTable->getNumOpacities(),
				     log( *temperatureIterator ),
				     log( targetDensity ), 
				     opacity );
		
		// The opacity vector contains the solution.  Now
		// we copy this solution into the OutputIterator
		for ( int i=0; i<ng; ++i, ++opacityIterator )
		    *opacityIterator = opacity[i];
		
	    }
	return opacityIterator;
    }

    /*!
     * \brief Opacity accessor that returns a vector of opacities that 
     *     corresponds to the provided temperature and density.
     */
Multigroup::OpacityType Multigroup::getOpacity(
    const double targetTemperature,
    const double targetDensity,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	int numGroups = spGandolfDataTable->getNumGroupBoundaries() - 1;
	OpacityType opacity( numGroups );
	// logorithmic interpolation:
	wrapper::wgintmglog( spGandolfDataTable->getLogTemperatures(),
			     spGandolfDataTable->getNumTemperatures(), 
			     spGandolfDataTable->getLogDensities(), 
			     spGandolfDataTable->getNumDensities(),
			     spGandolfDataTable->getNumGroupBoundaries(),
			     spGandolfDataTable->getLogOpacities(), 
			     spGandolfDataTable->getNumOpacities(),
			     log(targetTemperature),
			     log(targetDensity), 
			     opacity );
	return opacity;
    }

    /*!
     * \brief Opacity accessor that returns a vector of vectors of
     *     opacities that correspond to the provided vector of
     *     temperatures and a single density value.
     */
std::vector< Multigroup::OpacityType > Multigroup::getOpacity(
    const std::vector<double>& targetTemperature,
    const double targetDensity,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	int ng = spGandolfDataTable->getNumGroupBoundaries() - 1;
	std::vector< OpacityType > opacity( targetTemperature.size() );
	for ( int i=0; i<targetTemperature.size(); ++i )
	    {
		opacity[i].resize(ng);
		// logorithmic interpolation:
		wrapper::wgintmglog( spGandolfDataTable->getLogTemperatures(),
				     spGandolfDataTable->getNumTemperatures(), 
				     spGandolfDataTable->getLogDensities(), 
				     spGandolfDataTable->getNumDensities(),
				     spGandolfDataTable->getNumGroupBoundaries(),
				     spGandolfDataTable->getLogOpacities(), 
				     spGandolfDataTable->getNumOpacities(),
				     log(targetTemperature[i]),
				     log(targetDensity), 
				     opacity[i] );
	    }
		return opacity;
    }

    /*!
     * \brief Opacity accessor that returns a vector of vectors of
     *     opacities that correspond to the provided vector of
     *     densities and a single temperature value.
     */
std::vector< Multigroup::OpacityType > Multigroup::getOpacity(
    const double targetTemperature,
    const std::vector<double>& targetDensity,
    const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
    {
	int numGroups = spGandolfDataTable->getNumGroupBoundaries() - 1;
	std::vector< OpacityType > opacity( targetDensity.size() );
	for ( int i=0; i<targetDensity.size(); ++i )
	    {
		opacity[i].resize(numGroups);
		// logorithmic interpolation:
		wrapper::wgintmglog( spGandolfDataTable->getLogTemperatures(),
				     spGandolfDataTable->getNumTemperatures(), 
				     spGandolfDataTable->getLogDensities(), 
				     spGandolfDataTable->getNumDensities(),
				     spGandolfDataTable->getNumGroupBoundaries(),
				     spGandolfDataTable->getLogOpacities(), 
				     spGandolfDataTable->getNumOpacities(),
				     log(targetTemperature),
				     log(targetDensity[i]),  
				     opacity[i] );
	    }
	return opacity;
    }

    // This routine was removed because of confusion surrounding the
    // ordering of temperature, density and multigroup indicies in the 
    // returned opacity container.

// std::vector< Multigroup::OpacityType > Multigroup::getOpacity(
//     const std::vector<double>& targetTemperature,
//     const std::vector<double>& targetDensity,
//     const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const
//     {
// 	int ng = spGandolfDataTable->getNumGroupBoundaries() - 1;
// 	int nt = targetTemperature.size();
// 	int nd = targetDensity.size();
// 	std::vector< OpacityType > opacity( nt*nd );
// 	for ( int i=0; i<nt; ++i )
// 	    for ( int j=0; j<nd; ++j )
// 		{
// 		    opacity[i*nt+j].resize(ng);
// 		    // logorithmic interpolation:
// 		    wrapper::wgintmglog( spGandolfDataTable->getLogTemperatures(),
// 					 spGandolfDataTable->getNumTemperatures(), 
// 					 spGandolfDataTable->getLogDensities(), 
// 					 spGandolfDataTable->getNumDensities(),
// 					 spGandolfDataTable->getNumGroupBoundaries(),
// 					 spGandolfDataTable->getLogOpacities(), 
// 					 spGandolfDataTable->getNumOpacities(),
// 					 log(targetTemperature[i]),
// 					 log(targetDensity[j]), 
// 					 opacity[i*nd+j] );
// 		}
// 	return opacity;
//     }

} // end namespace rtt_cdi_gandolf


//---------------------------------------------------------------------------//
//                              end of GandolfOpacity.cc
//---------------------------------------------------------------------------//
