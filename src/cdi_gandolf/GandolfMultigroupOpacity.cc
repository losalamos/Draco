//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfMultigroupOpacity.cc
 * \author Kelly Thompson
 * \date   Mon Jan 22 15:24210 2001
 * \brief  GandolfMultigroupOpacity templated class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfMultigroupOpacity.hh"

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
     * \brief Constructor for GandolfMultigroupOpacity object.
     * 
     * See GandolfMultigroupOpacity.hh for details.
     */
    GandolfMultigroupOpacity::GandolfMultigroupOpacity( 
	const rtt_dsxx::SP< const GandolfFile >& in_spGandolfFile,
	int in_materialID,
	rtt_cdi::Model in_opacityModel,
	rtt_cdi::Reaction in_opacityReaction )
	: spGandolfFile( in_spGandolfFile ),
	materialID( in_materialID ),
	numKeys( 0 ),
	opacityModel( in_opacityModel ),
	opacityReaction( in_opacityReaction ),
	energyPolicyDescriptor( "mg" )
	{
	    // Verify that the requested material ID is available in the
	    // specified IPCRESS file.
	    if ( ! spGandolfFile->materialFound( materialID ) )
		throw gkeysException( -1 );
	    
	    // Retrieve keys available fo this material from the IPCRESS
	    // file.  wgkeys() returns vKnownKeys, numKeys and errorCode.
	    int errorCode = 
		wrapper::wgkeys( spGandolfFile->getDataFilename(),
				 materialID, vKnownKeys,
				 wrapper::maxKeys, numKeys );
	    if ( errorCode !=0 ) throw gkeysException( errorCode );
	    
	    // Create the data table object and fill it with the table
	    // data from the IPCRESS file.
	    spGandolfDataTable = new GandolfDataTable(
		energyPolicyDescriptor, opacityModel, opacityReaction,
		vKnownKeys, materialID, spGandolfFile );
	    
	} // end of GandolfData constructor
    
    /*!
     * \brief Default GandolfOpacity() destructor.
     *
     * \sa This is required to correctly release memory when a 
     *     GandolfMultigroupOpacity is destroyed.  This constructor's
     *     * definition must be declared in the implementation file so
     *     that * we can avoid including too many header files
     */
    GandolfMultigroupOpacity::~GandolfMultigroupOpacity()
	{
	    // empty
	}
    
    
    // --------- //
    // Accessors //
    // --------- //
    
    /*!
     * \brief Returns a "plain English" description of the opacity
     *     data that this class references. (e.g. "Multigroup Rosseland
     *     Scattering".) 
     *
     * The definition of this function is not included here to prevent 
     *     the inclusion of the GandolfFile.hh definitions within this 
     *     header file.
     */
    std::string GandolfMultigroupOpacity::getDataDescriptor() const 
	{
	    // call the correct function from the GandolfDataTable
	    // object. 
	    return spGandolfDataTable->getDataDescriptor(); 
	}
    
    /*!
     * \brief Returns the name of the associated IPCRESS file.
     *
     * The definition of this function is not included here to prevent 
     *     the inclusion of the GandolfFile.hh definitions within this 
     *     header file.
     */
    std::string GandolfMultigroupOpacity::getDataFilename() const 
	{
	    return spGandolfFile->getDataFilename(); 
	}
    
    /*!
     * \brief Opacity accessor that returns a single opacity (or a
     *     vector of opacities for the multigroup EnergyPolicy) that 
     *     corresponds to the provided temperature and density.
     */
    std::vector< double > GandolfMultigroupOpacity::getOpacity(
	double targetTemperature,
	double targetDensity ) const
	{ 
	    // number of groups in this multigroup set.
	    const int numGroups = spGandolfDataTable->getNumGroupBoundaries() - 1;
	    
	    // temporary opacity vector used by the wrapper.  The returned 
	    // data will be copied into the opacityIterator.
	    std::vector<double> opacity( numGroups );
	    
	    // logorithmic interpolation:
	    opacity = wrapper::wgintmglog( 
		spGandolfDataTable->getLogTemperatures(),
		spGandolfDataTable->getNumTemperatures(), 
		spGandolfDataTable->getLogDensities(), 
		spGandolfDataTable->getNumDensities(),
		spGandolfDataTable->getNumGroupBoundaries(),
		spGandolfDataTable->getLogOpacities(), 
		spGandolfDataTable->getNumOpacities(),
		log( targetTemperature ),
		log( targetDensity ) );
	    return opacity;
	}
    
    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     temperatures and a single density value.
     */
    std::vector< std::vector< double > > GandolfMultigroupOpacity::getOpacity(
	const std::vector<double>& targetTemperature,
	double targetDensity ) const
	{ 
	    int numGroups = spGandolfDataTable->getNumGroupBoundaries() - 1;
	    std::vector< std::vector< double > > opacity( targetTemperature.size() );
	    for ( int i=0; i<targetTemperature.size(); ++i )
		{
		    opacity[i].resize( numGroups );
		    // logorithmic interpolation:
		    opacity[i] = wrapper::wgintmglog( 
			spGandolfDataTable->getLogTemperatures(),
			spGandolfDataTable->getNumTemperatures(), 
			spGandolfDataTable->getLogDensities(), 
			spGandolfDataTable->getNumDensities(),
			spGandolfDataTable->getNumGroupBoundaries(),
			spGandolfDataTable->getLogOpacities(), 
			spGandolfDataTable->getNumOpacities(),
			log( targetTemperature[i] ),
			log( targetDensity ) );
		}
	    return opacity;
	}
    
    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     densities and a single temperature value.
     */
    std::vector< std::vector< double > > GandolfMultigroupOpacity::getOpacity(
	double targetTemperature,
	const std::vector<double>& targetDensity ) const
	{ 
	    int numGroups = spGandolfDataTable->getNumGroupBoundaries() - 1;
	    std::vector< std::vector< double > > opacity( targetDensity.size() );
	    for ( int i=0; i<targetDensity.size(); ++i )
		{
		    opacity[i].resize(numGroups);
		    // logorithmic interpolation:
		    opacity[i] = wrapper::wgintmglog( 
			spGandolfDataTable->getLogTemperatures(),
			spGandolfDataTable->getNumTemperatures(), 
			spGandolfDataTable->getLogDensities(), 
			spGandolfDataTable->getNumDensities(),
			spGandolfDataTable->getNumGroupBoundaries(),
			spGandolfDataTable->getLogOpacities(), 
			spGandolfDataTable->getNumOpacities(),
			log( targetTemperature ),
			log( targetDensity[i] ) );
		}
	    return opacity;
	}
    
    /*!
     * \brief Returns a vector of temperatures that define the cached
     *     opacity data table.
     */
    std::vector< double > GandolfMultigroupOpacity::getTemperatureGrid() const
	{
	    return spGandolfDataTable->getTemperatures();
	}
    
    /*!
     * \brief Returns the size of the temperature grid.
     */
    int GandolfMultigroupOpacity::getNumTemperatures() const
	{
	    return spGandolfDataTable->getNumTemperatures();
	}
    
    /*!
     * \brief Returns a vector of densities that define the cached
     *     opacity data table.
     */
    std::vector<double> GandolfMultigroupOpacity::getDensityGrid() const
	{
	    return spGandolfDataTable->getDensities();
	}
    
    /*! 
     * \brief Returns the size of the density grid.
     */
    int GandolfMultigroupOpacity::getNumDensities() const
	{
	    return spGandolfDataTable->getNumDensities();
	}
    
    /*!
     * \brief Returns a vector of energy values (keV) that define the
     *     energy boundaries of the cached multigroup opacity data
     *     table.  (This accessor is only valid for the Multigroup
     *     EnergyPolicy version of GandolfMultigroupOpacity.)
     */
    std::vector< double > GandolfMultigroupOpacity::getGroupBoundaries() const
	{
	    return spGandolfDataTable->getGroupBoundaries();
	}
    
    /*!
     * \brief Returns the number of group boundaries found in the
     *     current multigroup data set.
     */
    int GandolfMultigroupOpacity::getNumGroupBoundaries() const
	{
	    return spGandolfDataTable->getNumGroupBoundaries();
	}
    
} // end namespace rtt_cdi_gandolf

//---------------------------------------------------------------------------//
// end of GandolfMultigroupOpacity.cc
//---------------------------------------------------------------------------//
