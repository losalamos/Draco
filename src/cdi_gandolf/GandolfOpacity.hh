//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfOpacity.hh
 * \author Kelly Thompson
 * \date   Wed Jul 12 16:11:55 2000
 * \brief  GandolfOpacity class header file (derived from cdi/Opacity)
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfOpacity_hh__
#define __cdi_gandolf_GandolfOpacity_hh__

#include <vector>
#include <string>

#include "ds++/SP.hh"

namespace rtt_cdi_gandolf
{
    // -------------------- //
    // Forward declarations //
    // -------------------- //

    class GandolfFile;
    class GandolfDataTable;

    // --------------------- //
    // Enumerated data types //
    // --------------------- //

    /*!
     * \brief Physics model used to compute the opacity values.
     */
    enum Model
    {
	Rosseland,
	Plank
    };

    /*!
     * \brief Opacity reaction type stored in this opacity object.
     */
    enum Reaction
    {
	Total,      /*!< Total opacity value (scattering plus absorption). */
	Absorption, /*!< Absorption cross sections only. */
	Scattering  /*!< Scattering cross sections only. */
    };

//===========================================================================//
/*!
 * \class GandolfOpacity
 *
 * \brief GandolfOpacity allows the client code to retrieve opacity
 *        data for a particular material.  Each GandolfOpacity object
 *        represents a specific type of data defined by 5 attributes:
 *        an IPCRESS File (via a GandolfFile object), a material
 *        identifier, an energy policy, a physics model and a reaction
 *        type.
 *
 * \sa  This is a concrete class derived from cdi/Opacity.  This
 *      class allows to client to access the data in IPCRESS files
 *      via the Gandolf libraries.
 *      <p>
 * This class is designed to be used in conjuction with the CDI.
 * The client code will create a GandolfOpacity object and use this
 * object as an argument during the CDI instantiation.  The purpose of 
 * this class is to provide a mechanism for accessing data in IPCRESS
 * files and works by calling the Gandolf library provided by X-5.
 * The GandolfOpacity constructor expects two arguments: the IPCRESS
 * file name and a material identifier.  Once constructed this object
 * allows the client to access any data found in the IPCRESS file for
 * that one material.  The client code will need to create a separate
 * GandolfOpacity object for each material that it need information
 * about. 
 * <p>
 * Additionally, this is class uses policy templates.  GandolfOpacity
 * is templated on an EnergyPolicy (Gray or Multigroup - see the
 * bottom of this file).  These templated classes are aliased
 * (typedef) as GrayOpacity and MultigroupOpacity so the client code
 * can declare a GandolfOpacity<EnergyPolicy> object by using this
 * aliases. (see test/tGandolfOpacity.cc for examples).
 * <p>
 * When instantiated, the GandolfOpacity<EnergyPolicy> object creates
 * a GandolfDataTable object.  The IPCRESS data is cached in this
 * table object.  When the client requests an opacity value at a
 * specified temperature and density the GandolfOpcity object calls
 * GANDOLF library routines interpolate on the data cached in the
 * GandolfDataTable object.  
 * <p>
 * Every templated form of GandolfOpacity<EnergyPolicy> must be
 * instantiated in GandolfOpacity_pt.cc.
 */

/*!
 * \example cdi_gandolf/test/tGandolfOpacity.cc
 *
 * Example of GandolfOpacity usage independent of CDI.  In this
 * example we construct a GandolfOpacity object for the material
 * Aluminum (matID=10001 in our example IPCRESS file).  We then use
 * the GandolfOpacity object to compute a Rosseland Gray opacity
 * value for a specified material temperature and density.  In a
 * similar fashion we also request the GandolfOpacity object to return 
 * a vector of multigroup opacities for a specified temperature and
 * density. 
 *
 */
//===========================================================================//

// Using Policy Classes.  EnergyPolicy is either Gray or Multigroup.  

template < class EnergyPolicy >
class GandolfOpacity : public EnergyPolicy // rtt_cdi::Opacity
{

    // NESTED CLASSES AND TYPEDEFS

    /*!
     * \brief Opacity type is defined in the EnergyPolicy template.
     *        It is "double" for Gray and "vector<double>" for
     *        Multigroup. 
     */
    typedef typename EnergyPolicy::OpacityType OpacityType;


    // DATA

    // ----------------------- //
    // Specify unique material //
    // ----------------------- //

    /*!
     * \brief DS++ Smart Pointer to a GandolfFile object.
     *     spGandolfFile acts as a hook to link this object to an
     *     IPCRESS file.
     */
    const rtt_dsxx::SP<GandolfFile> spGandolfFile;

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
    std::vector<std::string> vKnownKeys;

    // --------------- //
    // Data specifiers //
    // --------------- //

    // the data table is specified by the EnergyPolicy (Gray or
    // Multigroup), the Model (Plank or Rosseland) and the Reaction
    // (Totak, Absorption, or Scattering).

    /*!
     * \brief The physics model that the current data set is based on.
     *        { Rosseland, Plank }.  This enumeration is defined
     *        above. 
     */
    const Model opacityModel;

    /*!
     * \brief The type of reaction rates that the current data set
     *        represents { Total, Scattering, Absorption }. This
     *        enumeration is defined above.
     */
    const Reaction opacityReaction;

    // -------------------- //
    // Opacity lookup table //
    // -------------------- //

    /*!
     * \brief spGandolfDataTable contains a cached copy of the
     *        requested IPCRESS opacity lookup table.
     *
     * There is a one-to-one relationship between GandolfOpacity and
     * GandolfDataTable. 
     */
    rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable;

  public:

    // ------------ //
    // Constructors //
    // ------------ //

    /*!
     * \brief This is the default GandolfOpacity constructor.  It
     *     requires four arguments plus the energy policy to be
     *     instantiated. 
     * 
     * \sa The combiniation of a data file and a material ID uniquely 
     *     specifies a material.  If we add the Model, Reaction and
     *     EnergyPolicy the opacity table is uniquely defined.
     *
     * \parameter _spGandolfFile This smart pointer links an IPCRESS
     *     file (via the GandolfFile object) to a GandolfOpacity
     *     object. There may be many GandolfOpacity objects per
     *     GandolfFile object but only one GandolfFile object for each 
     *     GandolfOpacity object.
     *
     * \parameter _materialID An identifier that links the
     *     GandolfOpacity object to a single material found in the
     *     specified IPCRESS file.
     *
     * \parameter _opacityModel The physics model that the current
     *     data set is based on.
     *
     * \parameter _opacityReaction The type of reaction rate that the
     *     current data set represents.
     */
    GandolfOpacity( const rtt_dsxx::SP<GandolfFile> _spGandolfFile,
		    const int _materialID, 
		    const Model _opacityModel,
		    const Reaction _opacityReaction );

    /*!
     * \brief Default GandolfOpacity() destructor.
     *
     * This is required to correctly release memory when a
     * GandolfOpacity<EnergyPolicy> is destroyed.
     */
    ~GandolfOpacity();

    // --------- //
    // Accessors //
    // --------- //

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
     * 
     * \parameter temperatureIterator The beginning position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter temperatureIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter densityIterator The beginning position of a STL
     *     container that holds a list of densities.
     *
     * \parameter densityIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     * 
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values corresponding to the
     *     given (temperature,density) tuple will be stored.
     *
     * \return A list (of type OutputIterator) of opacities are
     *     returned.  These opacities correspond to the temperature
     *     and density values provied in the two InputIterators.
     */
    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity( InputIterator temperatureIterator, 
			       InputIterator temperatureIteratorEnd,
			       InputIterator densityIterator, 
			       InputIterator densityIteratorEnd,
			       OutputIterator opacityIterator ) const;

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
     *
     * \parameter temperatureIterator The beginning position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter temperatureIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter targetDensity The single density value used when
     *     computing opacities for each given temperature.
     * 
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values (corresponding to the
     *     provided temperature and density values) will be stored.
     *
     * \return A list (of type OutputIterator) of opacities are
     *     returned.  These opacities correspond to the temperature
     *     and density values provied.
     */
    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity( InputIterator temperatureIterator,
			       InputIterator temperatureIteratorEnd,
			       const double targetDensity,
			       OutputIterator opacityIterator ) const;

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
     *
     * \parameter targetTemperature The single temperature value used when
     *     computing opacities for each given density.
     * 
     * \parameter densityIterator The beginning position of a STL
     *     container that holds a list of densities.
     *
     * \parameter densityIteratorEnd The end position of a STL
     *     container that holds a list of densities.
     *
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values (corresponding to the
     *     provided temperature and density values) will be stored.
     *
     * \return A list (of type OutputIterator) of opacities are
     *     returned.  These opacities correspond to the temperature
     *     and density values provied.
     */
    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity( const double targetTemperature,
			       InputIterator densityIterator, 
			       InputIterator densityIteratorEnd,
			       OutputIterator opacityIterator ) const;

    /*!
     * \brief Opacity accessor that returns a single opacity (or a
     *     vector of opacities for the multigroup EnergyPolicy) that 
     *     corresponds to the provided temperature and density.
     *
     * \parameter targetTemperature The temperature value for which an
     *     opacity value is being requested.
     *
     * \parameter targetDensity The density value for which an opacity 
     *     value is being requested.
     *
     * \return A single opacity (or a vector of opacities for the
     *     multigroup EnergyPolicy).
     */
    OpacityType getOpacity( const double targetTemperature,
			    const double targetDensity ) const; 

    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     temperatures and a single density value.
     *
     * \parameter targetTemperature A vector of temperature values for
     *     which opacity values are being requested.
     *
     * \parameter targetDensity The density value for which an opacity 
     *     value is being requested.
     *
     * \return A vector of opacities (or a vector of vectors of
     *     opacities for the multigroup EnergyPolicy).
     */
    std::vector< OpacityType > getOpacity( 
	const std::vector<double>& targetTemperature,
	const double targetDensity ) const; 

    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     densities and a single temperature value.
     *
     * \parameter targetTemperature The temperature value for which an 
     *     opacity value is being requested.
     *
     * \parameter targetDensity A vector of density values for which
     *     opacity values are being requested.
     *
     * \return A vector of opacities (or a vector of vectors of
     *     opacities for the multigroup EnergyPolicy).
     */
    std::vector< OpacityType > getOpacity( 
	const double targetTemperature,
	const std::vector<double>& targetDensity ) const; 

    // It is not clear how to assume order of opacity(temp,dens) when
    // accessed in this manner --> for now use the STL-style accessor
    // or a loop over one of the other vector-accessors.

//     std::vector< OpacityType > getOpacity( 
// 	const std::vector<double>& targetTemperature,
// 	const std::vector<double>& targetDensity ) const;

    /*!
     * \brief Returns a string that describes the templated
     *     EnergyPolicy.  Currently this will return either "mg" or
     *     "gray."
     */ 
    const std::string& getEnergyPolicyDescriptor() const { 
	return EnergyPolicy::getEnergyPolicyDescriptor(); };

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
    const std::string& getDataFilename() const {
	return spGandolfFile->getDataFilename(); };

    /*!
     * \brief Returns a vector of temperatures that define the cached
     *     opacity data table.
     */
    std::vector<double> getTemperatureGrid() const;

    /*!
     * \brief Returns a vector of densities that define the cached
     *     opacity data table.
     */
    std::vector<double> getDensityGrid() const;

    /*!
     * \brief Returns a vector of energy values (keV) that define the
     *     energy boundaries of the cached multigroup opacity data
     *     table.  (This accessor is only valid for the Multigroup
     *     EnergyPolicy version of GandolfOpacity.)
     */
    const std::vector<double>& getGroupBoundaries() const;
    
    /*!
     * \brief Returns the size of the temperature grid.
     */
    int getNumTemperatures() const;

    /*! 
     * \brief Returns the size of the density grid.
     */
    int getNumDensities() const;

    /*!
     * \brief Returns the number of group boundaries found in the
     *     current multigroup data set.
     */
    int getNumGroupBoundaries() const;

}; // end of class GandolfOpacity



//===========================================================================
//===========================================================================
//===========================================================================

// Should Multigroup and Gray classes be moved to the CDI package?

//===========================================================================
//===========================================================================
//===========================================================================





//===========================================================================//
/*!
 * \class Multigroup
 *
 * \brief This is an EnergyPolicy Class.  GnadolfOpacity will be
 *        derived from this class if it needs to manipulate multigroup
 *        data. 
 *
 * \sa The Multigroup EnergyPolicy sets the OpacityType to be
 * "vector<double>" 
 */
//===========================================================================//

class Multigroup
{
    // DATA

    /*!
     * \brief A string that describes the energy Policy ("mg").
     */
    const std::string energyPolicyDescriptor;

  protected:

    // -------- //
    // Typedefs //
    // -------- //
    
    /*!
     * \brief For multigroup data the Opacity type will be a vector of 
     *     doubles.  That is each "opacity" will be considered to be
     *     the set of multigroup opacity values corresponding to a
     *     temperature and density.
     *
     * This value must be private (not public) so that
     *     GandolfOpacity<EP> can see it. 
     */   
    typedef std::vector<double> OpacityType;

    // ---------------------------- //
    // Constructors and Destructors //
    // ---------------------------- //

    /*!
     * \brief This is the default Multigroup constructor.  It
     *     initializes the energyPolicyDescriptor to specify
     *     multigroup data.
     */
    Multigroup() 
	: energyPolicyDescriptor( "mg" ) { };

    /*!
     * \brief Default Multigroup descructor.  This is required to
     *     correctly release memory when this object is destroyed.
     */
    ~Multigroup() {};

    // --------- //
    // Accessors //
    // --------- //

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
     * 
     * \parameter temperatureIterator The beginning position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter temperatureIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter densityIterator The beginning position of a STL
     *     container that holds a list of densities.
     *
     * \parameter densityIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     * 
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values corresponding to the
     *     given (temperature,density) tuple will be stored.  The
     *     length of this container should be equal to the number of
     *     temperatures times the number of energy groups.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A list (of type OutputIterator) of opacities are
     *     returned.  These opacities correspond to the temperature
     *     and density values provied in the two InputIterators.
     */
    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity(
	InputIterator temperatureIterator, 
	InputIterator temperatureIteratorEnd,
	InputIterator densityIterator, 
	InputIterator densityIteratorEnd,
	OutputIterator opacityIterator,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable )
	const;

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of temperatures in an STL container.
     *     An set of opacity values will be returned for each temperature
     *     provided.  The opacity iterator should have length equal to 
     *     the length of the temperatureIterator times the number of
     *     energy groups.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     *
     * \parameter temperatureIterator The beginning position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter temperatureIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter targetDensity The single density value used when
     *     computing opacities for each given temperature.
     * 
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values (corresponding to the
     *     provided temperature and density values) will be stored.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A list (of type OutputIterator) of opacities is
     *     returned.  These opacities correspond to the temperature
     *     and density values provied.
     */
    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity(
	InputIterator temperatureIterator, 
	InputIterator temperatureIteratorEnd,
	const double targetDensity,
	OutputIterator opacityIterator,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const;

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of densities in an STL container.
     *     An opacity value will be returned for each density
     *     provided.  The opacity iterator should have length equal to 
     *     to the length of the densityIterator times the number of
     *     energy groups.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     *
     * \parameter targetTemperature The single temperature value used when
     *     computing opacities for each given density.
     * 
     * \parameter densityIterator The beginning position of a STL
     *     container that holds a list of densities.
     *
     * \parameter densityIteratorEnd The end position of a STL
     *     container that holds a list of densities.
     *
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values (corresponding to the
     *     provided temperature and density values) will be stored.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A list (of type OutputIterator) of opacities is
     *     returned.  These opacities correspond to the temperature
     *     and density values provied.
     */
    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity(
	const double targetTemperature,
	InputIterator densityIterator, 
	InputIterator densityIteratorEnd,
	OutputIterator opacityIterator,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const;

    /*!
     * \brief Opacity accessor that returns a vector of opacities that 
     *     corresponds to the provided temperature and density.
     *
     * \parameter targetTemperature The temperature value for which an
     *     opacity value is being requested.
     *
     * \parameter targetDensity The density value for which an opacity 
     *     value is being requested.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A multigroup opacities (a vector of group averaged
     *     opacities). 
     */
    OpacityType getOpacity( 
	const double targetTemperature,
	const double targetDensity,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable )
	const;

    /*!
     * \brief Opacity accessor that returns a vector of vectors of
     *     opacities that correspond to the provided vector of
     *     temperatures and a single density value.
     *
     * \parameter targetTemperature A vector of temperature values for
     *     which opacity values are being requested.
     *
     * \parameter targetDensity The density value for which an opacity 
     *     value is being requested.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A vector of multigroup opacities (each entry is a
     *     vector over energy groups). 
     */
    std::vector< OpacityType > getOpacity( 
	const std::vector<double>& targetTemperature,
	const double targetDensity,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable )
	const;

    /*!
     * \brief Opacity accessor that returns a vector of vectors of
     *     opacities that correspond to the provided vector of
     *     densities and a single temperature value.
     *
     * \parameter targetTemperature The temperature value for which an 
     *     opacity value is being requested.
     *
     * \parameter targetDensity A vector of density values for which
     *     opacity values are being requested.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A vector opacities (each entry is a vector over energy
     *     groups).
     */
    std::vector< OpacityType > getOpacity( 
	const double targetTemperature,
	const std::vector<double>& targetDensity,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable )
	const;

    // This routine was removed because of confusion surrounding the
    // ordering of temperature, density and multigroup indicies in the 
    // returned opacity container.

//     std::vector< OpacityType > getOpacity( 
// 	const std::vector<double>& targetTemperature,
// 	const std::vector<double>& targetDensity,
// 	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const;

    /*!
     * \brief Returns a string that describes the templated
     *     EnergyPolicy.  Currently this will return either "mg" or
     *     "gray."
     */ 
    const std::string& getEnergyPolicyDescriptor() const {
	return energyPolicyDescriptor; };

}; // end of class Multigroup


//===========================================================================//
/*!
 * \class Gray
 *
 * \brief This is an EnergyPolicy Class.  GnadolfOpacity will be
 *        derived from this class if it needs to manipulate multigroup
 *        data. 
 *
 * \sa The Gray EnergyPolicy sets the OpacityType to be "double".
 */
//===========================================================================//

class Gray
{
    // DATA

    /*!
     * \brief A string that describes the energy Policy ("gray").
     */
    const std::string energyPolicyDescriptor;

  protected:

    // -------- //
    // Typedefs //
    // -------- //
    
    /*!
     * \brief For gray data the Opacity type will be a double.  That
     *     is each "opacity" will be considered to be a single gray
     *     opacity value corresponding to a temperature and density.
     *
     * This value must be protected (not private) so that
     *     GandolfOpacity<EP> can see it. 
     */ 
    typedef double OpacityType;

    // ---------------------------- //
    // Constructors and Destructors //
    // ---------------------------- //

    /*!
     * \brief This is the default gray constructor.  It
     *     initializes the energyPolicyDescriptor to specify
     *     gray data.
     */
    Gray()
	: energyPolicyDescriptor( "gray" ) { };

    /*!
     * \brief Default Multigroup descructor.  This is required to
     *     correctly release memory when this object is destroyed.
     */
    ~Gray() { };

    // --------- //
    // Accessors //
    // --------- //

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
     * 
     * \parameter temperatureIterator The beginning position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter temperatureIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter densityIterator The beginning position of a STL
     *     container that holds a list of densities.
     *
     * \parameter densityIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     * 
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values corresponding to the
     *     given (temperature,density) tuple will be stored.  The
     *     length of this container should be equal to the number of
     *     temperatures provided by the temperatureIterator.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A list (of type OutputIterator) of opacities are
     *     returned.  These opacities correspond to the temperature
     *     and density values provied in the two InputIterators.
     */
    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity(
	InputIterator temperatureIterator, 
	InputIterator temperatureIteratorEnd,
	InputIterator densityIterator, 
	InputIterator densityIteratorEnd,
	OutputIterator opacityIterator,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable )
	const;

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of temperatures in an STL container.
     *     An set of opacity values will be returned for each temperature
     *     provided.  The opacity iterator should have length equal to 
     *     the length of the temperatureIterator.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     *
     * \parameter temperatureIterator The beginning position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter temperatureIteratorEnd The end position of a STL
     *     container that holds a list of temperatures.
     *
     * \parameter targetDensity The single density value used when
     *     computing opacities for each given temperature.
     * 
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values (corresponding to the
     *     provided temperature and density values) will be stored.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A list (of type OutputIterator) of opacities is
     *     returned.  These opacities correspond to the temperature
     *     and density values provied.
     */
    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity(
	InputIterator temperatureIterator, 
	InputIterator temperatureIteratorEnd,
	const double targetDensity,
	OutputIterator opacityIterator,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const;

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of densities in an STL container.
     *     An opacity value will be returned for each density
     *     provided.  The opacity iterator should have length equal to 
     *     to the length of the densityIterator.
     *
     * The InputIterator and OutputIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in GandolfOpacity_pt.cc.
     *
     * \parameter targetTemperature The single temperature value used when
     *     computing opacities for each given density.
     * 
     * \parameter densityIterator The beginning position of a STL
     *     container that holds a list of densities.
     *
     * \parameter densityIteratorEnd The end position of a STL
     *     container that holds a list of densities.
     *
     * \parameter opacityIterator The beginning position of a STL
     *     container into which opacity values (corresponding to the
     *     provided temperature and density values) will be stored.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A list (of type OutputIterator) of opacities is
     *     returned.  These opacities correspond to the temperature
     *     and density values provied.
     */
    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity(
	const double targetTemperature,
	InputIterator densityIterator, 
	InputIterator densityIteratorEnd,
	OutputIterator opacityIterator,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable )
	const;

    /*!
     * \brief Opacity accessor that returns a single opacity that 
     *     corresponds to the provided temperature and density.
     *
     * \parameter targetTemperature The temperature value for which an
     *     opacity value is being requested.
     *
     * \parameter targetDensity The density value for which an opacity 
     *     value is being requested.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A single gray opacity.
     */
    OpacityType getOpacity( 
	const double targetTemperature,
	const double targetDensity,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const;

    /*!
     * \brief Opacity accessor that returns a vector of
     *     opacities that correspond to the provided vector of
     *     temperatures and a single density value.
     *
     * \parameter targetTemperature A vector of temperature values for
     *     which opacity values are being requested.
     *
     * \parameter targetDensity The density value for which an opacity 
     *     value is being requested.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A vector of gray opacities.
     */
    std::vector< OpacityType > getOpacity( 
	const std::vector<double>& targetTemperature,
	const double targetDensity,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const;

    /*!
     * \brief Opacity accessor that returns a vector of
     *     opacities that correspond to the provided vector of
     *     densities and a single temperature value.
     *
     * \parameter targetTemperature The temperature value for which an 
     *     opacity value is being requested.
     *
     * \parameter targetDensity A vector of density values for which
     *     opacity values are being requested.
     *
     * \parameter spGandolfDataTable This object is associated with
     *     GandolfOpacity<EP> and contains a cached copy of the
     *     IPCRESS opacity lookup table.
     *
     * \return A vector of gray opacities.
     */
    std::vector< OpacityType > getOpacity( 
	const double targetTemperature,
	const std::vector<double>& targetDensity,
	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable )
	const;

    // This routine was removed because of confusion surrounding the
    // ordering of temperature, density and multigroup indicies in the 
    // returned opacity container.

//     std::vector< OpacityType > getOpacity( 
// 	const std::vector<double>& targetTemperature,
// 	const std::vector<double>& targetDensity,
// 	const rtt_dsxx::SP<GandolfDataTable> spGandolfDataTable ) const;

    /*!
     * \brief Returns a string that describes the templated
     *     EnergyPolicy.  Currently this will return either "mg" or
     *     "gray."
     */ 
    const std::string& getEnergyPolicyDescriptor() const {
	return energyPolicyDescriptor; };
};



//===========================================================================//
//                          MORE TYPEDEFS                      
//===========================================================================//

/*!
 * \brief An alias for GandolfOpacity<Gray>
 */
typedef GandolfOpacity<Gray> GrayOpacity;

/*!
 * \brief An alias for GandolfOpacity<Multigroup>
 */
typedef GandolfOpacity<Multigroup> MGOpacity;

/*!
 * \brief An alias for GandolfOpacity<Multigroup>
 */
typedef GandolfOpacity<Multigroup> MultigroupOpacity;

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfOpacity_hh__

//---------------------------------------------------------------------------//
//                end of cdi_gandolf/GandolfOpacity.hh
//---------------------------------------------------------------------------//
