//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_ipcress/IpcressDataTable.hh
 * \author Kelly Thompson
 * \date   Thu Oct 12 09:39:22 2000
 * \brief  Header file for IpcressDataTable
 * \note   Copyright (C) 2001-2010 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_ipcress_IpcressDataTable_hh__
#define __cdi_ipcress_IpcressDataTable_hh__

// We must include IpcressOpacity.hh so that the following two
// enumerated items are defined: { Model, Reaction } 
#include "cdi/OpacityCommon.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

namespace rtt_cdi_ipcress
{
    
// forward declaration (we don't need to include the full
// definition of IpcressFile).
    
class IpcressFile;
    
//===========================================================================//
/*!
 * \class IpcressDataTable encapsulates all of the data associated
 * with a specific opacity type (e.g.: total, plank, multigroup) for a 
 * single single material.
 *
 * When the user instantiates a IpcressOpacity object a
 * IpcressDataTable object is also created.  There is a one-to-one
 * correspondence between these two objects.  The 
 * IpcressDataTable object will load a single opacity table from the
 * IPCRESS file specified by the associated IpcressFile object.  The
 * table that is loaded is specified by the combination of 
 * { opacityModel, opacityReaction and the opacityEnergyDescriptor }.
 */
//===========================================================================//

class IpcressDataTable 
{

    // NESTED CLASSES AND TYPEDEFS

    // Data Descriptors:

    /*!
     * \brief A string that specifies the type of data being stored.  Possible
     *     values are rgray, ragray, rsgray, etc.  This key is provided to the
     *     Ipcress libraries as a data specifier.
     */
    std::string mutable ipcressDataTypeKey;

    /*!
     * \brief A string that specifies the type of data being stored.  This
     *     variables holds an English version of ipcressDataTypeKey.
     */
    std::string mutable dataDescriptor;

    /*!
     * \brief A string that specifies the energy model for the data being
     *     stored.  Possible values are "mg" or "gray".
     */
    std::string const opacityEnergyDescriptor;

    /*!
     * \brief An enumerated value defined in IpcressOpacity.hh that specifies
     *     the data model.  Possible values are "Rosseland" or "Plank".
     */
    rtt_cdi::Model const opacityModel;

    /*!
     * \brief An enumerated valued defined in IpcressOpacity.hh that specifies
     *     the reaction model.  Possible values are "Total", "Absorption" or
     *     "Scattering".
     */
    rtt_cdi::Reaction const opacityReaction;

    /*!
     * \brief A list of keys that are known by the IPCRESS file.
     */
    std::vector< std::string > const & vKnownKeys;
    
    /*!
     * \brief The IPCRESS material number assocated with the data
     *     contained in this object.
     */
    size_t const matID;
    
    /*!
     * \brief The IpcressFile object assocaiated with this data.
     */
    rtt_dsxx::SP< const IpcressFile > const spIpcressFile;


    // Data Sizes:

    // This object holds a single data table.  This table is loaded
    // from the IPCRSS file (specified by the IpcressFile object).
    // The size of this table is stored in the following variables.

    //! The number of temperature columns in the opacity table.
    size_t mutable numTemperatures;

    //! The number of density columns in the opacity table.
    size_t mutable numDensities;

    /*!
     * \brief The number of energy group boundaries in the opacity
     *     table (this entry is not used for gray data).
     */
    size_t mutable numGroupBoundaries;

    /*
     * \brief The number of entries in the opacity table.  This should 
     *     be equal to numTemperatures * numDensities *
     *     (numGroupBoundaries - 1).
     */
    size_t mutable numOpacities;


    // Data Tables:

    //! The temperature grid for this data set.
    std::vector<double> mutable logTemperatures;
    std::vector<double> mutable temperatures;

    //! The density grid for this data set.
    std::vector<double> mutable logDensities;
    std::vector<double> mutable densities;

    //! The energy group boundary grid for this data set.
    std::vector<double> mutable groupBoundaries;

    //! The opacity data table.
    std::vector<double> mutable logOpacities;

  public:

    // CREATORS

    /*!
     * \brief Standard IpcressDataTable constructor.
     *
     * The constructor requires that the data state to be
     *     completely defined.  With this information the DataTypeKey
     *     is set, then the data table sizes are loaded and finally
     *     the table data is loaded.
     *
     * \param opacityEnergyDescriptor This string variable 
     *     specifies the energy model { "gray" or "mg" } for the
     *     opacity data contained in this IpcressDataTable object. 
     * \param opacityModel This enumerated value specifies the
     *     physics model { Rosseland or Plank } for the opacity data
     *     contained in this object.  The enumeration is defined in
     *     IpcressOpacity.hh 
     * \param opacityReaction This enumerated value specifies the 
     *     interaction model { total, scattering, absorption " for the 
     *     opacity data contained in this object.  The enumeration is
     *     defined in IpcressOpacity.hh
     * \param vKnownKeys This vector of strings is a list of
     *     data keys that the IPCRESS file knows about.  This list is
     *     read from the IPCRESS file when a IpcressOpacity object is
     *     instantiated but before the associated IpcressDataTable
     *     object is created. 
     * \param matID The material identifier that specifies a
     *     particular material in the IPCRESS file to associate with
     *     the IpcressDataTable container.
     * \param spIpcressFile A DS++ SmartPointer to a IpcressFile
     *     object.  One GanolfFile object should exist for each
     *     IPCRESS file.  Many IpcressOpacity (and thus
     *     IpcressDataTable) objects may point to the same IpcressFile 
     *     object. 
     */    
    IpcressDataTable( std::string                const & opacityEnergyDescriptor,
		      rtt_cdi::Model                     opacityModel, 
		      rtt_cdi::Reaction                  opacityReaction,
		      std::vector< std::string > const & vKnownKeys,
		      size_t                             matID,
		      rtt_dsxx::SP< const IpcressFile > const & spIpcressFile );

    // ACCESSORS

    //! Retrieve the size of the temperature grid.
    size_t getNumTemperatures() const { return numTemperatures; };

    //! Retrieve the size of the density grid.
    size_t getNumDensities() const { return numDensities; };

    //! Retrieve the size of the energy boundary grid.
    size_t getNumGroupBoundaries() const { return numGroupBoundaries; };

    //! Retrieve the size of the opacity grid.
    size_t getNumOpacities() const { return numOpacities; };

    //! Retrieve the logarithmic temperature grid.
    const std::vector<double>& getLogTemperatures() const { 
	return logTemperatures; };
    const std::vector<double>& getTemperatures() const { 
	return temperatures; };

    //! Retrieve the logarithmic density grid.
    const std::vector<double>& getLogDensities() const {
	return logDensities; };
    const std::vector<double>& getDensities() const {
	return densities; };

    //! Retrieve the logarithmic opacity grid.
    const std::vector<double>& getLogOpacities() const {
	return logOpacities; };

    //! Retrieve the energy boundary grid.
    std::vector<double> const & getGroupBoundaries() const {
	return groupBoundaries; };

    //! Return a "plain English" description of the data table.
    std::string const & getDataDescriptor() const {
	return dataDescriptor; };

    /*!
     * \brief Return a "plain English" description of the energy
     *     policy. 
     */
    std::string const & getEnergyPolicyDescriptor() const {
	return opacityEnergyDescriptor; };

    //! Perform linear interploation of log(opacity) values.
    double interpOpac( double const T, double const rho ) const;
    
  private:

    /*!
     * \brief This function sets both "ipcressDataTypeKey" and
     *     "dataDescriptor" based on the values given for
     *     opacityEnergyDescriptor, opacityModel and opacityReaction.
     */
    void setIpcressDataTypeKey() const;

    /*!
     * \brief Load the table sizes from the IPCRESS file and resize
     *     the vector containers for the actual data tables. 
     */
    void setIpcressDataTableSizes() const;

    /*!
     * \brief Load the temperature, density, energy boundary and
     *     opacity opacity tables from the IPCRESS file.  Convert all
     *     tables (except energy boundaries) to log values.
     */
    void loadDataTable();

    /*!
     * \brief Search "keys" for "key".  If found return true,
     *     otherwise return false.
     */
    template < typename T >
    bool key_available( const T &key, const std::vector<T> &keys ) const;

};

} // end namespace rtt_cdi_ipcress

#endif // __cdi_ipcress_IpcressDataTable_hh__

//---------------------------------------------------------------------------//
// end of cdi_ipcress/IpcressDataTable.hh
//---------------------------------------------------------------------------//
