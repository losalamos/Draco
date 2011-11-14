//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_ipcress/IpcressDataTable.cc
 * \author Kelly Thompson
 * \date   Thu Oct 12 09:39:22 2000
 * \brief  Implementation file for IpcressDataTable objects.
 * \note   Copyright (C) 2001-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "IpcressDataTable.hh"  // the associated header file.
#include "IpcressFile.hh"       // we have a SP to a GandofFile object.
#include "cdi/OpacityCommon.hh" // defines Model and Reaction
                                // enumerated values.
#include "ds++/Assert.hh"
#include <cmath>                // we need to define log(double)


// ------------------------- //
// NAMESPACE RTT_CDI_IPCRESS //
// ------------------------- //

namespace rtt_cdi_ipcress
{

/*!
 * \brief IpcressData Table constructor.
 *
 *     The constructor requires that the data state be completely defined.
 *     With this information the DataTypeKey is set, then the data table sizes
 *     are loaded and finally the table data is loaded.
 *
 * \param opacityEnergyDescriptor This string variable specifies the energy
 *     model { "gray" or "mg" } for the opacity data contained in this
 *     IpcressDataTable object.
 * \param opacityModel This enumerated value specifies the physics model {
 *     Rosseland or Planck } for the opacity data contained in this object.
 *     The enumeration is defined in IpcressOpacity.hh
 * \param opacityReaction This enumerated value specifies the interaction
 *     model { total, scattering, absorption " for the opacity data contained
 *     in this object.  The enumeration is defined in IpcressOpacity.hh
 * \param vKnownKeys This vector of strings is a list of data keys that the
 *     IPCRESS file knows about.  This list is read from the IPCRESS file when
 *     a IpcressOpacity object is instantiated but before the associated
 *     IpcressDataTable object is created.
 * \param matID The material identifier that specifies a particular material
 *     in the IPCRESS file to associate with the IpcressDataTable container.
 * \param spIpcressFile A ds++ SmartPointer to a IpcressFile object.  One
 *     GanolfFile object should exist for each IPCRESS file.  Many
 *     IpcressOpacity (and thus IpcressDataTable) objects may point to the
 *     same IpcressFile object.
 */
IpcressDataTable::IpcressDataTable( 
    std::string              const & in_opacityEnergyDescriptor,
    rtt_cdi::Model                   in_opacityModel, 
    rtt_cdi::Reaction                in_opacityReaction,
    std::vector<std::string> const & in_vKnownKeys,
    int                              in_matID,
    rtt_dsxx::SP< const IpcressFile > const & in_spIpcressFile )
    : ipcressDataTypeKey( "" ),
      dataDescriptor( "" ),
      opacityEnergyDescriptor ( in_opacityEnergyDescriptor ),
      opacityModel( in_opacityModel ),
      opacityReaction( in_opacityReaction ),
      vKnownKeys ( in_vKnownKeys ),
      matID ( in_matID ),
      spIpcressFile( in_spIpcressFile ),
      numTemperatures( 0 ),
      numDensities( 0 ),
      numGroupBoundaries( 0 ),
      numOpacities( 0 ),
      logTemperatures(),
      temperatures(),
      logDensities(),
      densities(),
      groupBoundaries(),
      logOpacities()
{
    // Obtain the Ipcress keyword for the opacity data type specified by the
    // EnergyPolicy, opacityModel and the opacityReaction.  Valid keywords
    // are: { ramg, rsmg, rtmg, pmg, rgray, ragray, rsgray, pgray } This
    // function also ensures that the requested data type is available in the
    // IPCRESS file.
    setIpcressDataTypeKey();
    
    // Retrieve the size of the data set and resize the vector containers.
    setIpcressDataTableSizes();
    
    // Retrieve table data (temperatures, densities, group boundaries and
    // opacities.  These are stored as logorithmic values.
    loadDataTable();
    
} // end of IpcressDataTable constructor.

    
// ----------------- //
// PRIVATE FUNCTIONS //
// ----------------- //


/*!
 * \brief This function sets both "ipcressDataTypeKey" and
 *     "dataDescriptor" based on the values given for
 *     opacityEnergyDescriptor, opacityModel and opacityReaction.
 */
void IpcressDataTable::setIpcressDataTypeKey( ) const
{
    // Build the Ipcress key for the requested data.  Valid
    // keys are:
    // { ramg, rsmg, rtmg, pmg, rgray, ragray, rsgray, pgray }

    if ( opacityEnergyDescriptor == "gray" )
    {
        switch ( opacityModel ) {
            case ( rtt_cdi::ROSSELAND ) :

                switch ( opacityReaction ) {
                    case ( rtt_cdi::TOTAL ) :
                        ipcressDataTypeKey = "rgray";
                        dataDescriptor = "Gray Rosseland Total";
                        break;
                    case ( rtt_cdi::ABSORPTION ) :
                        ipcressDataTypeKey = "ragray";
                        dataDescriptor = "Gray Rosseland Absorption";
                        break;
                    case ( rtt_cdi::SCATTERING ) :
                        // *** NOTE: THIS KEY DOES NOT ACTUALLY EVER EXIST *** //
                        // See LA-UR-01-5543
                        ipcressDataTypeKey = "rsgray";
                        dataDescriptor = "Gray Rosseland Scattering";
                        break;
                    default :
                        Assert(false);
                        break;
                }
                break;

            case ( rtt_cdi::PLANCK ) :
			
                switch ( opacityReaction ) {
                    case ( rtt_cdi::TOTAL ) :
                        // *** NOTE: THIS KEY DOES NOT ACTUALLY EVER EXIST *** //
                        // See LA-UR-01-5543
                        ipcressDataTypeKey = "ptgray";
                        dataDescriptor = "Gray Planck Total";
                        break;
                    case ( rtt_cdi::ABSORPTION ) :
                        ipcressDataTypeKey = "pgray";
                        dataDescriptor = "Gray Planck Absorption";
                        break;
                    case ( rtt_cdi::SCATTERING ) :
                        // *** NOTE: THIS KEY DOES NOT ACTUALLY EVER EXIST *** //
                        // See LA-UR-01-5543
                        ipcressDataTypeKey = "psgray";
                        dataDescriptor = "Gray Planck Scattering";
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
            case ( rtt_cdi::ROSSELAND ) :

                switch ( opacityReaction ) {
                    case ( rtt_cdi::TOTAL ) :
                        ipcressDataTypeKey = "rtmg";
                        dataDescriptor = "Multigroup Rosseland Total";
                        break;
                    case ( rtt_cdi::ABSORPTION ) :
                        ipcressDataTypeKey = "ramg";
                        dataDescriptor = "Multigroup Rosseland Absorption";
                        break;
                    case ( rtt_cdi::SCATTERING ) :
                        ipcressDataTypeKey = "rsmg";
                        dataDescriptor = "Multigroup Rosseland Scattering";
                        break;
                    default :
                        Assert(false);
                        break;
                }
                break;
			
            case ( rtt_cdi::PLANCK ) :
			
                switch ( opacityReaction ) {
                    case ( rtt_cdi::TOTAL ) :
                        // *** NOTE: THIS KEY DOES NOT ACTUALLY EVER EXIST *** //
                        // See LA-UR-01-5543
                        ipcressDataTypeKey = "ptmg";
                        dataDescriptor = "Multigroup Planck Total";
                        break;
                    case ( rtt_cdi::ABSORPTION ) :
                        ipcressDataTypeKey = "pmg";
                        dataDescriptor = "Multigroup Planck Absorption";
                        break;
                    case ( rtt_cdi::SCATTERING ) :
                        // *** NOTE: THIS KEY DOES NOT ACTUALLY EVER EXIST *** //
                        // See LA-UR-01-5543
                        ipcressDataTypeKey = "psmg";
                        dataDescriptor = "Multigroup Planck Scattering";
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
    Insist( key_available( ipcressDataTypeKey, vKnownKeys ),
            "requested opacity type is not available in the IPCRESS file.");

}

/*!
 * \brief Load the temperature, density, energy boundary and
 *     opacity opacity tables from the IPCRESS file.  Convert all
 *     tables (except energy boundaries) to log values.
 */
void IpcressDataTable::setIpcressDataTableSizes() const
{
    // int idum = 0;
    int errorCode = 0;
    // A different wrapper routine must be called for
    // multigroup and gray data.  We choose the correct
    // wrapper by comparing the opacityEnergyDescriptor.
    if ( opacityEnergyDescriptor == "mg" )
    {
        // Returns: numTemperatures, numDensities,
        // numGroupBoundaries and numOpacities.
        errorCode = 0;
            // wrapper::wgchgrids(
            //     spIpcressFile->getDataFilename(),
            //     matID, numTemperatures, numDensities,
            //     numGroupBoundaries, idum, numOpacities );
        Check(false);
    }
    else // gray
    {
        // Returns: numTemperatures, numDensities,
        // numGroupBoundaries and numOpacities.
        errorCode = 0;
            // wrapper::wgchgrids(
            //     spIpcressFile->getDataFilename(),
            //     matID, numTemperatures, numDensities,
            //     numGroupBoundaries, numOpacities, idum );
        Check(false);
    }
	    
    // if the wrapper returned an error code the we need to
    // throw an exception.
    // if ( errorCode != 0 ) throw gchgridsException( errorCode );

    // Resize the data containers based on the newly loaded
    // size parameters.
    temperatures.resize( numTemperatures );
    logTemperatures.resize( numTemperatures );
    densities.resize( numDensities );
    logDensities.resize( numDensities );
    groupBoundaries.resize( numGroupBoundaries );
    logOpacities.resize( numOpacities );	    
}

/*!
 * \brief Load the temperature, density, energy boundary and
 *     opacity opacity tables from the IPCRESS file.  Convert all
 *     tables (except energy boundaries) to log values.
 */
void IpcressDataTable::loadDataTable() const
{
    int errorCode = 0;
    // A different wrapper routine must be called for
    // multigroup and gray data.  We choose the correct
    // wrapper by comparing the opacityEnergyDescriptor.
    if ( opacityEnergyDescriptor == "mg" )
    {
        // Returns: logTemperatures, logDensities,
        // groupBoundaries and logOpacities.
        errorCode = 0;
            // wrapper::wggetmg( 
            //     spIpcressFile->getDataFilename(), 
            //     matID, ipcressDataTypeKey,
            //     temperatures, numTemperatures,
            //     densities, numDensities,
            //     groupBoundaries, numGroupBoundaries,
            //     logOpacities, numOpacities );
        // if ( errorCode != 0 )
        //throw ggetmgException( errorCode );
                Check(false);

    }
    else // "gray"
    {
        // Returns: logTemperatures, logDensities and
        // logOpacities.
        errorCode = 0;
            // wrapper::wggetgray( 
            //     spIpcressFile->getDataFilename(), 
            //     matID, ipcressDataTypeKey,
            //     temperatures, numTemperatures,
            //     densities, numDensities,
            //     logOpacities, numOpacities );
        // if the wrapper returned an error code the we need to
        // throw an exception.
        //if ( errorCode != 0 )
        //throw ggetgrayException( errorCode );
        Check(false);
    }
    // The interpolation routines expect everything to be in
    // log form so we only store the logorithmic temperature,
    // density and opacity data.
    for ( int i=0; i<numTemperatures; ++i )
        logTemperatures[i] = std::log( temperatures[i] );
    for ( int i=0; i<numDensities; ++i )
        logDensities[i] = std::log( densities[i] );
    for ( int i=0; i<numOpacities; ++i )
        logOpacities[i] = std::log( logOpacities[i] );
}

/*! 
 * \brief This function returns "true" if "key" is found in the list
 *        of "keys".  This is a static member function.
 */
template < typename T >
bool IpcressDataTable::key_available( const T &key, 
                                      const std::vector<T> &keys ) const
{
    // Loop over all available keys.  If the requested key
    // matches one in the list return true.  If we reach the end
    // of the list without a match return false.
    for ( size_t i=0; i<keys.size(); ++i )
        if ( key == keys[i] ) return true;
    return false;
	    
} // end of IpcressDataTable::key_available( string, vector<string> )
    
} // end namespace rtt_cdi_ipcress

//---------------------------------------------------------------------------//
// end of IpcressDataTable.cc
//---------------------------------------------------------------------------//
