//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfGrayOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 22 14:11:10 2001
 * \brief  GandolfGrayOpacity templated class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfGrayOpacity.hh"

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
 * \brief Constructor for GandolfGrayOpacity object.
 * 
 * See GandolfGrayOpacity.hh for details.
 */
GandolfGrayOpacity::GandolfGrayOpacity( 
    const rtt_dsxx::SP< GandolfFile > _spGandolfFile,
    const int _materialID,
    const rtt_cdi::Model _opacityModel,
    const rtt_cdi::Reaction _opacityReaction )
    : spGandolfFile( _spGandolfFile ),
      materialID( _materialID ),
      numKeys( 0 ),
      opacityModel( _opacityModel ),
      opacityReaction( _opacityReaction ),
      energyPolicyDescriptor( "gray" )
    {
 	// Verify that the requested material ID is available in the
 	// specified IPCRESS file.
 	if ( ! spGandolfFile->materialFound( materialID ) )
 	    throw gkeysException( -1 );

	// Retrieve keys available for this material from the IPCRESS
	// file.  wgkeys() returns vKnownKeys, numKeys and errorCode.
	int errorCode;
	wrapper::wgkeys( spGandolfFile->getDataFilename(), materialID, 
			 vKnownKeys, wrapper::maxKeys, numKeys,
			 errorCode );
	if ( errorCode !=0 ) throw gkeysException( errorCode );

	// Create the data table object and fill it with the table
	// data from the IPCRESS file.
 	spGandolfDataTable = new GandolfDataTable(
 	    energyPolicyDescriptor,
	    opacityModel, 
	    opacityReaction,
	    vKnownKeys,
	    materialID, 
	    spGandolfFile );

    } // end of GandolfData constructor

/*!
 * \brief Desctructor for GandolfGrayOpacity class.
 */ 
GandolfGrayOpacity::~GandolfGrayOpacity()
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
const std::string& GandolfGrayOpacity::getDataDescriptor() const 
    {
	// call the correct function from the GandolfDataTable
	// object. 
	return spGandolfDataTable->getDataDescriptor(); 
    }

     /*!
     * \brief Returns the name of the associated IPCRESS file.
     *
     *     The definition of this function is not included here to
     *     prevent the inclusion of the GandolfFile.hh definitions
     *     within this header file.
     */
    const std::string& GandolfGrayOpacity::getDataFilename() const 
	{ 
	    return spGandolfFile->getDataFilename(); 
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
     *     been done for a few STL containers in GandolfGrayOpacity_pt.cc.
     */
// template < class EnergyPolicy >
// template < class InputIterator, class OutputIterator >
// OutputIterator 
// GandolfGrayOpacity< EnergyPolicy >::getOpacity(
//     InputIterator temperatureIterator, 
//     InputIterator temperatureIteratorEnd,
//     InputIterator densityIterator, 
//     InputIterator densityIteratorEnd,
//     OutputIterator opacityIterator ) const
//     { 
// 	// call the correct accessor from the PolicyTemplate
// 	return EnergyPolicy::getOpacity(
// 	    temperatureIterator, temperatureIteratorEnd,
// 	    densityIterator, densityIteratorEnd,
// 	    opacityIterator, spGandolfDataTable );
//     }

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
     *     been done for a few STL containers in GandolfGrayOpacity_pt.cc.
     */
// template < class EnergyPolicy >
// template < class InputIterator, class OutputIterator >
// OutputIterator 
// GandolfGrayOpacity< EnergyPolicy >::getOpacity(
//     InputIterator temperatureIterator, 
//     InputIterator temperatureIteratorEnd,
//     const double targetDensity,
//     OutputIterator opacityIterator ) const
//     { 
// 	// call the correct accessor from the PolicyTemplate
// 	return EnergyPolicy::getOpacity(
// 	    temperatureIterator, temperatureIteratorEnd,
// 	    targetDensity,
// 	    opacityIterator, spGandolfDataTable );
//     }

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
     *     been done for a few STL containers in GandolfGrayOpacity_pt.cc.
     */
// template < class EnergyPolicy >
// template < class InputIterator, class OutputIterator >
// OutputIterator 
// GandolfGrayOpacity< EnergyPolicy >::getOpacity(
//     const double targetTemperature,
//     InputIterator densityIterator, 
//     InputIterator densityIteratorEnd,
//     OutputIterator opacityIterator ) const
//     { 
// 	// call the correct accessor from the PolicyTemplate
// 	return EnergyPolicy::getOpacity(
// 	    targetTemperature,
// 	    densityIterator, densityIteratorEnd,
// 	    opacityIterator, spGandolfDataTable );
//     }

    /*!
     * \brief Opacity accessor that returns a single opacity (or a
     *     vector of opacities for the multigroup EnergyPolicy) that 
     *     corresponds to the provided temperature and density.
     */
double GandolfGrayOpacity::getOpacity(
    const double targetTemperature,
    const double targetDensity ) const
    { 
	double opacity;
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
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     temperatures and a single density value.
     */
std::vector< double > GandolfGrayOpacity::getOpacity(
    const std::vector<double>& targetTemperature,
    const double targetDensity ) const
    { 
	std::vector< double > opacity( targetTemperature.size() );
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
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     densities and a single temperature value.
     */
std::vector< double > GandolfGrayOpacity::getOpacity(
    const double targetTemperature,
    const std::vector<double>& targetDensity ) const
    { 
	std::vector< double > opacity( targetDensity.size() );
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

    /*!
     * \brief Returns a vector of temperatures that define the cached
     *     opacity data table.
     */
const std::vector< double >& GandolfGrayOpacity::getTemperatureGrid() const
    {
	return spGandolfDataTable->getTemperatures();
    }

    /*!
     * \brief Returns the size of the temperature grid.
     */
int GandolfGrayOpacity::getNumTemperatures() const
    {
	return spGandolfDataTable->getNumTemperatures();
    }

    /*!
     * \brief Returns a vector of densities that define the cached
     *     opacity data table.
     */
const std::vector<double>& GandolfGrayOpacity::getDensityGrid() const
    {
	return spGandolfDataTable->getDensities();
    }

    /*! 
     * \brief Returns the size of the density grid.
     */
int GandolfGrayOpacity::getNumDensities() const
    {
	return spGandolfDataTable->getNumDensities();
    }

} // end namespace rtt_cdi_gandolf

//---------------------------------------------------------------------------//
// end of GandolfGrayOpacity.cc
//---------------------------------------------------------------------------//
