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
	const rtt_dsxx::SP< const GandolfFile >& in_spGandolfFile,
	int in_materialID,
	rtt_cdi::Model in_opacityModel,
	rtt_cdi::Reaction in_opacityReaction )
	: spGandolfFile( in_spGandolfFile ),
	materialID( in_materialID ),
	numKeys( 0 ),
	opacityModel( in_opacityModel ),
	opacityReaction( in_opacityReaction ),
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
    std::string GandolfGrayOpacity::getDataDescriptor() const 
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
    std::string GandolfGrayOpacity::getDataFilename() const 
	{ 
	    return spGandolfFile->getDataFilename(); 
	}
 
    /*!
     * \brief Opacity accessor that returns a single opacity (or a
     *     vector of opacities for the multigroup EnergyPolicy) that 
     *     corresponds to the provided temperature and density.
     */
    double GandolfGrayOpacity::getOpacity(
	double targetTemperature,
	double targetDensity ) const
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
	double targetDensity ) const
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
	double targetTemperature,
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
    std::vector< double > GandolfGrayOpacity::getTemperatureGrid() const
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
    std::vector<double> GandolfGrayOpacity::getDensityGrid() const
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
