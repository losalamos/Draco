//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfDataTable.cc
 * \author Kelly Thompson
 * \date   Thu Oct 12 09:39:22 2000
 * \brief  Implementation file for GandolfDataTable objects.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "GandolfDataTable.hh"  // the associated header file.

#include "GandolfOpacity.hh"    // defines Model and Reaction
                                // enumerated values.

#include "GandolfWrapper.hh"    // we call the wrapper routines.

#include "GandolfException.hh"  // must catch gandolf exceptions when
                                // calling wrapper routines.

#include "GandolfFile.hh"       // we have a SP to a GandofFile object.

#include "ds++/Assert.hh"

#include <cmath>                // we need to define log(double)


// ------------------------- //
// NAMESPACE RTT_CDI_GANDOLF //
// ------------------------- //

namespace rtt_cdi_gandolf
{

    /*!
     * \brief GandolfData Table constructor.
     *
     * \sa The constructor requires that the data state be completely
     *     defined.  With this information the DataTypeKey is set,
     *     then the data table sizes are loaded and finally the table
     *     data is loaded. 
     *
     * \parameter _opacityEnergyDescriptor This string variable 
     *     specifies the energy model { "gray" or "mg" } for the
     *     opacity data contained in this GandolfDataTable object. 
     *     
     * \parameter _opacityModel This enumerated value specifies the
     *     physics model { Rosseland or Plank } for the opacity data
     *     contained in this object.  The enumeration is defined in
     *     GandolfOpacity.hh 
     *
     * \parameter _opacityReaction This enumerated value specifies the 
     *     interaction model { total, scattering, absorption " for the 
     *     opacity data contained in this object.  The enumeration is
     *     defined in GandolfOpacity.hh
     *
     * \parameter _vKnownKeys This vector of strings is a list of
     *     data keys that the IPCRESS file knows about.  This list is
     *     read from the IPCRESS file when a GandolfOpacity object is
     *     instantiated but before the associated GandolfDataTable
     *     object is created. 
     *
     * \parameter _matID The material identifier that specifies a
     *     particular material in the IPCRESS file to associate with
     *     the GandolfDataTable container.
     *
     * \parameter _spGandolfFile A DS++ SmartPointer to a GandolfFile
     *     object.  One GanolfFile object should exist for each
     *     IPCRESS file.  Many GandolfOpacity (and thus
     *     GandolfDataTable) objects may point to the same GandolfFile 
     *     object.  
     */
    GandolfDataTable::GandolfDataTable( 
 	const std::string _opacityEnergyDescriptor,
 	const Model _opacityModel, 
	const Reaction _opacityReaction,
	const std::vector<std::string>& _vKnownKeys,
	const int _matID,
	const rtt_dsxx::SP<GandolfFile> _spGandolfFile )
	: opacityEnergyDescriptor ( _opacityEnergyDescriptor ),
	  opacityModel( _opacityModel ),
	  opacityReaction( _opacityReaction ),
	  vKnownKeys ( _vKnownKeys ),
	  matID ( _matID ),
	  spGandolfFile( _spGandolfFile ),
  	  gandolfDataTypeKey( "" ),
	  dataDescriptor( "" ),
	  numTemperatures( 0 ),
	  numDensities( 0 ),
	  numGroupBoundaries( 0 ),
	  numOpacities( 0 )
	{
	    // Obtain the Gandolf keyword for the opacity data type
	    // specified by the EnergyPolicy, opacityModel and the
	    // opacityReaction.  Valid keywords are:
	    // { ramg, rsmg, rtmg, pmg, rgray, ragray, rsgray, pgray }
	    // This function also ensures that the requested data type 
	    // is available in the IPCRESS file.
	    setGandolfDataTypeKey();

	    // Retrieve the size of the data set and resize the vector 
	    // containers.
	    setGandolfDataTableSizes();

	    // Retrieve table data (temperatures, densities, group
	    // boundaries and opacities.  These are stored as
	    // logorithmic values.
	    loadDataTable();

	} // end of GandolfDataTable constructor.

    
// ----------------- //
// PRIVATE FUNCTIONS //
// ----------------- //


    /*!
     * \brief This function sets both "gandolfDataTypeKey" and
     *     "dataDescriptor" based on the values given for
     *     opacityEnergyDescriptor, opacityModel and opacityReaction.
     */
    void GandolfDataTable::setGandolfDataTypeKey( )
	{
	    // Build the Gandolf key for the requested data.  Valid
	    // keys are:
	    // { ramg, rsmg, rtmg, pmg, rgray, ragray, rsgray, pgray }

	    if ( opacityEnergyDescriptor == "gray" )
		{
		    switch ( opacityModel ) {
		    case ( Rosseland ) :

			switch ( opacityReaction ) {
			case ( Total ) :
			    gandolfDataTypeKey = "rgray";
			    dataDescriptor = "Gray Rosseland Total";
			    break;
			case ( Absorption ) :
			    gandolfDataTypeKey = "ragray";
			    dataDescriptor = "Gray Rosseland Absorption";
			    break;
			case ( Scattering ) :
			    gandolfDataTypeKey = "rsgray";
			    dataDescriptor = "Gray Rosseland Scattering";
			    break;
			default :
			    Assert(false);
			    break;
			}
			break;

		    case ( Plank ) :
			
			switch ( opacityReaction ) {
			case ( Total ) :
			    gandolfDataTypeKey = "pgray";
			    dataDescriptor = "Gray Plank Total";
			    break;
			case ( Absorption ) :
			    gandolfDataTypeKey = "pagray";
			    dataDescriptor = "Gray Plank Absorption";
			    break;
			case ( Scattering ) :
			    gandolfDataTypeKey = "psgray";
			    dataDescriptor = "Gray Plank Scattering";
			    break;
			default :
			    Assert(false);
			    break;
			}
			break;
			
		    default :
			Assert(false);
			break;
		    }
		}
	    else // "mg"
		{
		    switch ( opacityModel ) {
		    case ( Rosseland ) :

			switch ( opacityReaction ) {
			case ( Total ) :
			    gandolfDataTypeKey = "rtmg";
			    dataDescriptor = "Multigroup Rosseland Total";
			    break;
			case ( Absorption ) :
			    gandolfDataTypeKey = "ramg";
			    dataDescriptor = "Multigroup Rosseland Absorption";
			    break;
			case ( Scattering ) :
			    gandolfDataTypeKey = "rsmg";
			    dataDescriptor = "Multigroup Rosseland Scattering";
			    break;
			default :
			    Assert(false);
			    break;
			}
			break;
			
		    case ( Plank ) :
			
			switch ( opacityReaction ) {
			case ( Total ) :
			    gandolfDataTypeKey = "pmg";
			    dataDescriptor = "Multigroup Plank Total";
			    break;
			case ( Absorption ) :
			    gandolfDataTypeKey = "pamg";
			    dataDescriptor = "Multigroup Plank Absorption";
			    break;
			case ( Scattering ) :
			    gandolfDataTypeKey = "psmg";
			    dataDescriptor = "Multigroup Plank Scattering";
			    break;
			default :
			    Assert(false);
			    break;
			}
			break;
			
		    default :
			Assert(false);
			break;
		    }
		}

	    // Verify that the requested opacity type is available in
	    // the IPCRESS file.
 	    if ( ! key_available( gandolfDataTypeKey, vKnownKeys ) )
		throw gkeysException( -2 );
	}

    /*!
     * \brief Load the temperature, density, energy boundary and
     *     opacity opacity tables from the IPCRESS file.  Convert all
     *     tables (except energy boundaries) to log values.
     */
    void GandolfDataTable::setGandolfDataTableSizes()
	{
	    int idum, errorCode;
	    // A different wrapper routine must be called for
	    // multigroup and gray data.  We choose the correct
	    // wrapper by comparing the opacityEnergyDescriptor.
	    if ( opacityEnergyDescriptor == "mg" )
		{
		    // Returns: numTemperatures, numDensities,
		    // numGroupBoundaries and numOpacities.
		    wrapper::wgchgrids(
			spGandolfFile->getDataFilename(),
			matID, numTemperatures, numDensities,
			numGroupBoundaries, idum, numOpacities,
			errorCode );
		}
	    else // gray
		{
		    // Returns: numTemperatures, numDensities,
		    // numGroupBoundaries and numOpacities.
		    wrapper::wgchgrids(
			spGandolfFile->getDataFilename(),
			matID, numTemperatures, numDensities,
			numGroupBoundaries, numOpacities, idum, 
			errorCode );
		}

	    // if the wrapper returned an error code the we need to
	    // throw an exception.
	    if ( errorCode != 0 ) throw gchgridsException( errorCode );

	    // Resize the data containers based on the newly loaded
	    // size parameters.
	    logTemperatures.resize( numTemperatures );
	    logDensities.resize( numDensities );
	    groupBoundaries.resize( numGroupBoundaries );
	    logOpacities.resize( numOpacities );	    
	}

    /*!
     * \brief Load the temperature, density, energy boundary and
     *     opacity opacity tables from the IPCRESS file.  Convert all
     *     tables (except energy boundaries) to log values.
     */
    void GandolfDataTable::loadDataTable()
	{
	    int errorCode;
	    // A different wrapper routine must be called for
	    // multigroup and gray data.  We choose the correct
	    // wrapper by comparing the opacityEnergyDescriptor.
	    if ( opacityEnergyDescriptor == "mg" )
		{
		    // Returns: logTemperatures, logDensities,
		    // groupBoundaries and logOpacities.
		    wrapper::wggetmg( 
			spGandolfFile->getDataFilename(), matID, gandolfDataTypeKey,
			logTemperatures, numTemperatures, numTemperatures,
			logDensities, numDensities, numDensities,
			groupBoundaries, numGroupBoundaries, numGroupBoundaries,
			logOpacities, numOpacities, numOpacities,
			errorCode );
		    if ( errorCode != 0 )
			throw ggetmgException( errorCode );
		}
	    else // "gray"
		{
		    // Returns: logTemperatures, logDensities and
		    // logOpacities.
		    wrapper::wggetgray( 
			spGandolfFile->getDataFilename(), matID, gandolfDataTypeKey,
			logTemperatures, numTemperatures, numTemperatures,
			logDensities, numDensities, numDensities,
			logOpacities, numOpacities, numOpacities,
			errorCode );
		    // if the wrapper returned an error code the we need to
		    // throw an exception.
		    if ( errorCode != 0 )
			throw ggetgrayException( errorCode );
		}
	    // The interpolation routines expect everything to be in
	    // log form so we only store the logorithmic temperature,
	    // density and opacity data.
	    for ( int i=0; i<numTemperatures; ++i )
		logTemperatures[i] = log( logTemperatures[i] );
	    for ( int i=0; i<numDensities; ++i )
		logDensities[i] = log( logDensities[i] );
	    for ( int i=0; i<numOpacities; ++i )
		logOpacities[i] = log( logOpacities[i] );
	}

/*! 
  * \brief This function returns "true" if "key" is found in the list
  *        of "keys".  This is a static member function.
  */
template < typename T >
bool GandolfDataTable::key_available( const T &key, 
				      const std::vector<T> &keys ) const
    {
	// Loop over all available keys.  If the requested key
	// matches one in the list return true.  If we reach the end
	// of the list without a match return false.
	for ( int i=0; i<keys.size(); ++i )
	    if ( key == keys[i] ) return true;
	return false;
	
    } // end of GandolfDataTable::key_available( string, vector<string> )

} // end namespace rtt_cdi_gandolf

//---------------------------------------------------------------------------//
// end of GandolfDataTable.cc
//---------------------------------------------------------------------------//
