//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfGrayOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 22 13:23:37 2001
 * \brief  GandolfGrayOpacity class header file (derived from cdi/GrayOpacity)
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfGrayOpacity_hh__
#define __cdi_gandolf_GandolfGrayOpacity_hh__

#include <vector>
#include <string>

#include "ds++/SP.hh"
#include "cdi/GrayOpacity.hh"

namespace rtt_cdi_gandolf
{
    // -------------------- //
    // Forward declarations //
    // -------------------- //

    class GandolfFile;
    class GandolfDataTable;

//===========================================================================//
/*!
 * \class GandolfGrayOpacity
 *
 * \brief provides access to gray opacity data located in IPCRESS files.
 *
 *        GandolfGrayOpacity allows the client code to retrieve opacity
 *        data for a particular material.  Each GandolfOpacity object
 *        represents a specific type of data defined by five
 *        attributes: an IPCRESS File (via a GandolfFile object), a
 *        material identifier, an energy model (already selected
 *        since this is a Gray Opacity class), a physics model and a
 *        reaction type.
 *
 *      This is a concrete class derived from cdi/GrayOpacity.  This
 *      class allows to client to access the data in IPCRESS files
 *      via the Gandolf libraries.
 * <p>
 *      This class is designed to be used in conjuction with the CDI.
 *      The client code will create a GandolfGrayOpacity object and
 *      use this object as an argument during the CDI instantiation.
 *      The purpose of this class is to provide a mechanism for
 *      accessing data in IPCRESS files and works by calling the
 *      Gandolf library provided by X-5.  The GandolfGrayOpacity
 *      constructor expects four arguments: a hook to IPCRESS data
 *      file (spGandolfFile), a material identifier, an opacity model
 *      (Rosseland or Plank) and an opacity reaction specifier (total,
 *      scattering or absorption).  Once constructed, this object
 *      allows the client to access any data found in the IPCRESS file
 *      for that one material.  The client code will need to create a
 *      separate GandolfGrayOpacity object for each material that it
 *      needs information about. Multiple opacity objects can exist
 *      per IPCRESS file.
 * <p> 
 *      This class only provides access to gray opacity data.  If the
 *      user needs multigroup opacity IPCRESS data he/she should use
 *      the cdi_gandolf/GandolfMultigroupOpacity class.
 * <p>
 *      When instantiated, the GandolfGrayOpacity object creates a
 *      GandolfDataTable object.  The IPCRESS data is cached in this
 *      table object.  When the client requests an opacity value at a
 *      specified temperature and density the GandolfGrayOpcity object
 *      calls the appropriate GANDOLF library routine, which in turn,
 *      interpolates on the data cached in the GandolfDataTable
 *      object.
 * <p>
 *      When compiling DRACO with support for the IPCRESS file reader
 *      (via Gandolf) you must add the following option on the
 *      configure line:<br><br>
 *      <tt>    --with-gandolf-lib=${VENDORS}/gandolf/IRIX64/lib64</tt><br><br>
 *      Currently, the Gandolf library is only available on 64 bit
 *      IRIX architectures.
 */

/*!
 * \example cdi_gandolf/test/tGandolfOpacity.cc
 *
 * Example of GandolfGrayOpacity usage independent of CDI.  In
 *      this example we construct a GandolfGrayOpacity object for the
 *      material Aluminum (matID=10001 in our example IPCRESS file).
 *      We then use the GandolfGrayOpacity object to compute a
 *      Rosseland Gray opacity value for a specified material,
 *      temperature and density.  Other forms of the getOpacity()
 *      accessor are tested along with accessors that return
 *      information about the data set and the cached data table.
 *
 * \example cdi_gandolf/test/tGandolfWithCDI.cc
 * 
 * This example tests and demonstrates how to use the cdi_gandolf
 * package as a plug-in for the CDI class.
 */
//===========================================================================//

class GandolfGrayOpacity : public rtt_cdi::GrayOpacity
{

    // DATA

    // ----------------------- //
    // Specify unique material //
    // ----------------------- //

    /*!
     * \brief DS++ Smart Pointer to a GandolfFile object.
     *     spGandolfFile acts as a hook to link this object to an
     *     IPCRESS file.
     */
    const rtt_dsxx::SP< GandolfFile > spGandolfFile;

    /*!
     * \brief Identification number for one of the materials found in
     *     the IPCRESS file pointed to by spGandolfFile.
     */
    const int materialID;

    // -------------------- //
    // Available data types //
    // -------------------- //
    
    // The IPCRESS file only holds specific data for each of its materials.

    /*!
     * \brief Number of types of data found in the IPCRESS file.
     */
    int numKeys;

    /*!
     * \brief A list of keys known by the IPCRESS file.
     */
    std::vector< std::string > vKnownKeys;

    // --------------- //
    // Data specifiers //
    // --------------- //

    /*!
     * \brief The physics model that the current data set is based on.
     *        { Rosseland, Plank }.  This enumeration is defined
     *        in cdi/OpacityCommon.hh.
     */
    const rtt_cdi::Model opacityModel;

    /*!
     * \brief The type of reaction rates that the current data set
     *        represents { Total, Scattering, Absorption }. This
     *        enumeration is defined in cdi/OpacityCommon.hh.
     */
    const rtt_cdi::Reaction opacityReaction;

    /*!
     * \brief A string that identifies the energy model for this
     *     class.
     */
    const std::string energyPolicyDescriptor;

    // -------------------- //
    // Opacity lookup table //
    // -------------------- //

    /*!
     * \brief spGandolfDataTable contains a cached copy of the
     *        requested IPCRESS opacity lookup table.
     *
     * There is a one-to-one relationship between GandolfGrayOpacity and
     * GandolfDataTable. 
     */
    rtt_dsxx::SP< GandolfDataTable > spGandolfDataTable;

  public:

    // ------------ //
    // Constructors //
    // ------------ //

    /*!
     * \brief This is the default GandolfGrayOpacity constructor.  It
     *     requires four arguments plus the energy policy (this class)
     *     to be instantiated.
     * 
     *     The combiniation of a data file and a material ID uniquely 
     *     specifies a material.  If we add the Model, Reaction and
     *     EnergyPolicy the opacity table is uniquely defined.
     *
     * \param _spGandolfFile This smart pointer links an IPCRESS
     *     file (via the GandolfFile object) to a GandolfOpacity
     *     object. There may be many GandolfOpacity objects per
     *     GandolfFile object but only one GandolfFile object for each 
     *     GandolfOpacity object.
     * \param _materialID An identifier that links the
     *     GandolfOpacity object to a single material found in the
     *     specified IPCRESS file.
     * \param _opacityModel The physics model that the current
     *     data set is based on.
     * \param _opacityReaction The type of reaction rate that the
     *     current data set represents. 
     */
    GandolfGrayOpacity( const rtt_dsxx::SP< GandolfFile > _spGandolfFile,
			const int _materialID, 
			const rtt_cdi::Model _opacityModel,
			const rtt_cdi::Reaction _opacityReaction );

    /*!
     * \brief Default GandolfOpacity() destructor.
     *
     *     This is required to correctly release memory when a
     *     GandolfGrayOpacity is destroyed.  We define the destructor
     *     in the implementation file to avoid including the
     *     unnecessary header files.
     */
    ~GandolfGrayOpacity();

    // --------- //
    // Accessors //
    // --------- //

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of (temperature,density) tuples.
     *     An opacity value will be returned for each tuple.  The
     *     temperature and density iterators are required to
     *     be the same length.  The opacity iterator should also have
     *     this same length.
     * 
     * \param temperatureFirst The beginning position of a STL
     *     container that holds a list of temperatures.
     * \param temperatureLast The end position of a STL
     *     container that holds a list of temperatures.
     * \param densityFirst The beginning position of a STL
     *     container that holds a list of densities.
     * \param densityLast
     *     container that holds a list of temperatures.
     * \param opacityFirst The beginning position of a STL
     *     container into which opacity values corresponding to the
     *     given (temperature,density) tuple will be stored.
     * \return A list (of type OpacityIterator) of opacities are
     *     returned.  These opacities correspond to the temperature
     *     and density values provied in the two InputIterators.
     */
    template < class TemperatureIterator, class DensityIterator,
               class OpacityIterator >
    OpacityIterator getOpacity( TemperatureIterator temperatureFirst,
				TemperatureIterator temperatureLast,
				DensityIterator densityFirst, 
				DensityIterator densityLast,
				OpacityIterator opacityFirst ) const;

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of temperatures in an STL container.
     *     An opacity value will be returned for each temperature
     *     provided.  The opacity iterator should be the same length
     *     as the temperature STL container.
     *
     * \param temperatureFirst The beginning position of a STL
     *     container that holds a list of temperatures.
     * \param temperatureLast The end position of a STL
     *     container that holds a list of temperatures.
     * \param targetDensity The single density value used when
     *     computing opacities for each given temperature.
     * \param opacityFirst The beginning position of a STL
     *     container into which opacity values (corresponding to the
     *     provided temperature and density values) will be stored.
     * \return A list (of type OpacityIterator) of opacities are
     *     returned.  These opacities correspond to the temperature
     *     and density values provied.
     */
    template < class TemperatureIterator, class OpacityIterator >
    OpacityIterator getOpacity( TemperatureIterator temperatureFirst,
				TemperatureIterator temperatureLast,
				const double targetDensity,
				OpacityIterator opacityFirst ) const;

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of densities in an STL container.
     *     An opacity value will be returned for each density
     *     provided.  The opacity iterator should be the same length
     *     as the density STL container.
     *
     * \param targetTemperature The single temperature value used when
     *     computing opacities for each given density.
     * \param densityFirst The beginning position of a STL
     *     container that holds a list of densities.
     * \param densityLast The end position of a STL
     *     container that holds a list of densities.
     * \param opacityFirst beginning position of a STL
     *     container into which opacity values (corresponding to the
     *     provided temperature and density values) will be stored.
     * \return A list (of type OpacityIterator) of opacities are
     *     returned.  These opacities correspond to the temperature
     *     and density values provied.
     */
    template < class DensityIterator, class OpacityIterator >
    OpacityIterator getOpacity( const double targetTemperature,
				DensityIterator densityFirst, 
				DensityIterator densityLast,
				OpacityIterator opacityFirst ) const;
    
    /*!
     * \brief Opacity accessor that returns a single opacity that 
     *     corresponds to the provided temperature and density.
     *
     * \param targetTemperature The temperature value for which an
     *     opacity value is being requested.
     * \param targetDensity The density value for which an opacity 
     *     value is being requested.
     * \return A single opacity.
     */
    double getOpacity( const double targetTemperature,
		       const double targetDensity ) const; 

    /*!
     * \brief Opacity accessor that returns a vector of opacities 
     *     that correspond to the provided vector of
     *     temperatures and a single density value.
     *
     * \param targetTemperature A vector of temperature values for
     *     which opacity values are being requested.
     * \param targetDensity The density value for which an opacity 
     *     value is being requested.
     * \return A vector of opacities.
     */
    std::vector< double > getOpacity( 
	const std::vector<double>& targetTemperature,
	const double targetDensity ) const; 

    /*!
     * \brief Opacity accessor that returns a vector of opacities that
     *     correspond to the provided vector of 
     *     densities and a single temperature value.
     * \param targetTemperature The temperature value for which an 
     *     opacity value is being requested.
     * \param targetDensity A vector of density values for which
     *     opacity values are being requested.
     * \return A vector of opacities.
     */
    std::vector< double > getOpacity( 
	const double targetTemperature,
	const std::vector<double>& targetDensity ) const; 

    // It is not clear how to assume order of opacity(temp,dens) when
    // accessed in this manner --> for now use the STL-style accessor
    // or a loop over one of the other vector-accessors.

//     std::vector< double > getOpacity( 
// 	const std::vector<double>& targetTemperature,
// 	const std::vector<double>& targetDensity ) const;

    /*!
     * \brief Returns a string that describes the templated
     *     EnergyPolicy.  Currently this will return either "mg" or
     *     "gray."
     */ 
    const std::string& getEnergyPolicyDescriptor() const {
	return energyPolicyDescriptor; };

    /*!
     * \brief Returns a "plain English" description of the opacity
     *     data that this class references. (e.g. "Gray Rosseland
     *     Scattering".) 
     *
     * The definition of this function is not included here to prevent 
     *     the inclusion of the GandolfFile.hh definitions within this 
     *     header file.
     */
    const std::string& getDataDescriptor() const;

    /*!
     * \brief Returns the name of the associated IPCRESS file.
     *
     * The definition of this function is not included here to prevent 
     *     the inclusion of the GandolfFile.hh definitions within this 
     *     header file.
     */
    const std::string& getDataFilename() const;

    /*!
     * \brief Returns a vector of temperatures that define the cached
     *     opacity data table.
     * 
     * We do not return a const reference because this function
     * must construct this information from more fundamental tables.
     */
    const std::vector< double >& getTemperatureGrid() const;

    /*!
     * \brief Returns a vector of densities that define the cached
     *     opacity data table.
     * 
     * We do not return a const reference because this function
     * must construct this information from more fundamental tables.
     */
    const std::vector< double >& getDensityGrid() const;

    /*!
     * \brief Returns the size of the temperature grid.
     */
    int getNumTemperatures() const;

    /*! 
     * \brief Returns the size of the density grid.
     */
    int getNumDensities() const;

}; // end of class GandolfGrayOpacity

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfGrayOpacity_hh__

//---------------------------------------------------------------------------//
//                end of cdi_gandolf/GandolfGrayOpacity.hh
//---------------------------------------------------------------------------//
