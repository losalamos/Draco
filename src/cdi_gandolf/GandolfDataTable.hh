//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfDataTable.hh
 * \author Kelly Thompson
 * \date   Thu Oct 12 09:39:22 2000
 * \brief  Header file for GandolfDataTable
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfDataTable_hh__
#define __cdi_gandolf_GandolfDataTable_hh__

// We must include GandolfOpacity.hh so that the following two
// enumerated items are defined: { Model, Reaction } 
#include "cdi/OpacityCommon.hh"

#include "ds++/SP.hh"

#include <vector>
#include <string>

namespace rtt_cdi_gandolf
{
    
    // forward declaration (we don't need to include the full
    // definition of GandolfFile).
    
    class GandolfFile;
    
//===========================================================================//
/*!
 * \class GandolfDataTable encapsulates all of the data associated
 * with a specific opacity type (e.g.: total, plank, multigroup) for a 
 * single single material.
 *
 * When the user instantiates a GandolfOpacity object a
 * GandolfDataTable object is also created.  There is a one-to-one
 * correspondence between these two objects.  The 
 * GandolfDataTable object will load a single opacity table from the
 * IPCRESS file specified by the associated GandolfFile object.  The
 * table that is loaded is specified by the combination of 
 * { opacityModel, opacityReaction and the opacityEnergyDescriptor }.
 */
//===========================================================================//

class GandolfDataTable 
{

    // NESTED CLASSES AND TYPEDEFS

    // Data Descriptors:

    /*!
     * \brief A string that specifies the type of data being stored.
     *     Possible values are rgray, ragray, rsgray, etc.  This key
     *     is provided to the Gandolf libraries as a data specifier.
     */
    mutable std::string gandolfDataTypeKey;

    /*!
     * \brief A string that specifies the type of data being stored.
     *     This variables holds an English version of
     *     gandolfDataTypeKey. 
     */
    mutable std::string dataDescriptor;

    /*!
     * \brief A string that specifies the energy model for the data
     *     being stored.  Possible values are "mg" or "gray".
     */
    const std::string opacityEnergyDescriptor;

    /*!
     * \brief An enumerated value defined in GandolfOpacity.hh that
     *     specifies the data model.  Possible values are "Rosseland"
     *     or "Plank".
     */
    const rtt_cdi::Model opacityModel;

    /*!
     * \brief An enumerated valued defined in GandolfOpacity.hh that
     *     specifies the reaction model.  Possible values are "Total", 
     *     "Absorption" or "Scattering".
     */
    const rtt_cdi::Reaction opacityReaction;

    /*!
     * \brief A list of keys that are known by the IPCRESS file.
     */
    const std::vector< std::string >& vKnownKeys;
    
    /*!
     * \brief The IPCRESS material number assocated with the data
     *     contained in this object.
     */
    const int matID;
    
    /*!
     * \brief The GandolfFile object assocaiated with this data.
     */
    const rtt_dsxx::SP< const GandolfFile > spGandolfFile;


    // Data Sizes:

    // This object holds a single data table.  This table is loaded
    // from the IPCRSS file (specified by the GandolfFile object).
    // The size of this table is stored in the following variables.

    /*!
     * \brief The number of temperature columns in the opacity table.
     */
    mutable int numTemperatures;

    /*!
     * \brief The number of density columns in the opacity table.
     */
    mutable int numDensities;

    /*!
     * \brief The number of energy group boundaries in the opacity
     *     table (this entry is not used for gray data).
    */
    mutable int numGroupBoundaries;

    /*
     * \brief The number of entries in the opacity table.  This should 
     *     be equal to numTemperatures * numDensities *
     *     (numGroupBoundaries - 1).
     */
    mutable int numOpacities;


    // Data Tables:

    /*!
     * \brief The temperature grid for this data set.
     */
    mutable std::vector<double> logTemperatures;
    mutable std::vector<double> temperatures;

    /*!
     * \brief The density grid for this data set.
     */
    mutable std::vector<double> logDensities;
    mutable std::vector<double> densities;

    /*!
     * \brief The energy group boundary grid for this data set.
     */
    mutable std::vector<double> groupBoundaries;

    /*!
     * \brief The opacity data table.
     */
    mutable std::vector<double> logOpacities;

  public:

    // CREATORS

    /*!
     * \brief Standard GandolfDataTable constructor.
     *
     * The constructor requires that the data state to be
     *     completely defined.  With this information the DataTypeKey
     *     is set, then the data table sizes are loaded and finally
     *     the table data is loaded.
     *
     * \param opacityEnergyDescriptor This string variable 
     *     specifies the energy model { "gray" or "mg" } for the
     *     opacity data contained in this GandolfDataTable object. 
     * \param opacityModel This enumerated value specifies the
     *     physics model { Rosseland or Plank } for the opacity data
     *     contained in this object.  The enumeration is defined in
     *     GandolfOpacity.hh 
     * \param opacityReaction This enumerated value specifies the 
     *     interaction model { total, scattering, absorption " for the 
     *     opacity data contained in this object.  The enumeration is
     *     defined in GandolfOpacity.hh
     * \param vKnownKeys This vector of strings is a list of
     *     data keys that the IPCRESS file knows about.  This list is
     *     read from the IPCRESS file when a GandolfOpacity object is
     *     instantiated but before the associated GandolfDataTable
     *     object is created. 
     * \param matID The material identifier that specifies a
     *     particular material in the IPCRESS file to associate with
     *     the GandolfDataTable container.
     * \param spGandolfFile A DS++ SmartPointer to a GandolfFile
     *     object.  One GanolfFile object should exist for each
     *     IPCRESS file.  Many GandolfOpacity (and thus
     *     GandolfDataTable) objects may point to the same GandolfFile 
     *     object. 
     */    
    GandolfDataTable( const std::string& opacityEnergyDescriptor,
		      rtt_cdi::Model opacityModel, 
		      rtt_cdi::Reaction opacityReaction,
		      const std::vector< std::string >& vKnownKeys,
		      int matID,
		      const rtt_dsxx::SP< const GandolfFile >& spGandolfFile );

    // ACCESSORS

    /*!
     * \brief Retrieve the size of the temperature grid.
     */
    int getNumTemperatures() const { return numTemperatures; };

    /*!
     * \brief Retrieve the size of the density grid.
     */
    int getNumDensities() const { return numDensities; };

    /*!
     * \brief Retrieve the size of the energy boundary grid.
     */
    int getNumGroupBoundaries() const { return numGroupBoundaries; };

    /*!
     * \brief Retrieve the size of the opacity grid.
     */
    int getNumOpacities() const { return numOpacities; };

    /*!
     * \brief Retrieve the logarithmic temperature grid.
     */
    const std::vector<double>& getLogTemperatures() const { 
	return logTemperatures; };
    const std::vector<double>& getTemperatures() const { 
	return temperatures; };

    /*!
     * \brief Retrieve the logarithmic density grid.
     */
    const std::vector<double>& getLogDensities() const {
	return logDensities; };
    const std::vector<double>& getDensities() const {
	return densities; };

    /*!
     * \brief Retrieve the logarithmic opacity grid.
     */
    const std::vector<double>& getLogOpacities() const {
	return logOpacities; };

    /*!
     * \brief Retrieve the energy boundary grid.
     */
    const std::vector<double>& getGroupBoundaries() const {
	return groupBoundaries; };

    /*!
     * \brief Return a "plain English" description of the data table.
     */
    const std::string& getDataDescriptor() const {
	return dataDescriptor; };

    /*!
     * \brief Return a "plain English" description of the energy
     *     policy. 
     */
    const std::string& getEnergyPolicyDescriptor() const {
	return opacityEnergyDescriptor; };

  private:

    /*!
     * \brief This function sets both "gandolfDataTypeKey" and
     *     "dataDescriptor" based on the values given for
     *     opacityEnergyDescriptor, opacityModel and opacityReaction.
     */
    void setGandolfDataTypeKey() const;

    /*!
     * \brief Load the table sizes from the IPCRESS file and resize
     *     the vector containers for the actual data tables. 
     */
    void setGandolfDataTableSizes() const;

    /*!
     * \brief Load the temperature, density, energy boundary and
     *     opacity opacity tables from the IPCRESS file.  Convert all
     *     tables (except energy boundaries) to log values.
     */
    void loadDataTable() const;

    /*!
     * \brief Search "keys" for "key".  If found return true,
     *     otherwise return false.
     */
    template < typename T >
    bool key_available( const T &key, const std::vector<T> &keys ) const;

};

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfDataTable_hh__

//---------------------------------------------------------------------------//
// end of cdi_gandolf/GandolfDataTable.hh
//---------------------------------------------------------------------------//
