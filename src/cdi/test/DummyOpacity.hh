//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyOpacity.hh
 * \author Kelly Thompson
 * \date   Wed Jul 12 16:11:55 2000
 * \brief  DummyOpacity class header file (derived from cdi/Opacity)
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_DummyOpacity_hh__
#define __cdi_DummyOpacity_hh__

#include <vector>
#include <string>

#include "ds++/SP.hh"
#include "../Opacity.hh"

namespace rtt_dummy_opacity
{

template < class EnergyPolicy >
class DummyOpacity : rtt_cdi::Opacity< EnergyPolicy >
{

    // NESTED CLASSES AND TYPEDEFS

    /*!
     * \brief Opacity type is defined in the EnergyPolicy template.
     *        It is "double" for Gray and "vector<double>" for
     *        Multigroup. 
     */
    typedef typename EnergyPolicy::OpacityType OpacityType;

    // DATA
    
    const std::string dataFilename;
    const std::string dataDescriptor;

    std::vector<double> groupBoundaries;

  public:

    // ------------ //
    // Constructors //
    // ------------ //

    DummyOpacity( );
    ~DummyOpacity( );

    // --------- //
    // Accessors //
    // --------- //

    const std::string& getDataDescriptor() const;

    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity( InputIterator temperatureIterator, 
			       InputIterator temperatureIteratorEnd,
			       InputIterator densityIterator, 
			       InputIterator densityIteratorEnd,
			       OutputIterator opacityIterator ) const;

    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity( InputIterator temperatureIterator,
			       InputIterator temperatureIteratorEnd,
			       const double targetDensity,
			       OutputIterator opacityIterator ) const;

    template < class InputIterator, class OutputIterator >
    OutputIterator getOpacity( const double targetTemperature,
			       InputIterator densityIterator, 
			       InputIterator densityIteratorEnd,
			       OutputIterator opacityIterator ) const;

    OpacityType getOpacity( const double targetTemperature,
			    const double targetDensity ) const; 

    std::vector< OpacityType > getOpacity( 
	const std::vector<double>& targetTemperature,
	const double targetDensity ) const; 

    std::vector< OpacityType > getOpacity( 
	const double targetTemperature,
	const std::vector<double>& targetDensity ) const; 

    const std::string& getEnergyPolicyDescriptor() const { 
	return EnergyPolicy::getEnergyPolicyDescriptor(); };

    const std::string& getDataFilename() const {
	return dataFilename; };

    std::vector<double> getTemperatureGrid() const;

    std::vector<double> getDensityGrid() const;

    const std::vector<double>& getGroupBoundaries() const;
    
    int getNumTemperatures() const;

    int getNumDensities() const;

    int getNumGroupBoundaries() const;

}; // end of class GandolfOpacity





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
	OutputIterator opacityIterator ) const;

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
	OutputIterator opacityIterator ) const;

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
	OutputIterator opacityIterator ) const;

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
	const double targetDensity ) const;

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
	const double targetDensity ) const;

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
	const std::vector<double>& targetDensity ) const;

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
	OutputIterator opacityIterator ) const;

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
	OutputIterator opacityIterator ) const;

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
	OutputIterator opacityIterator ) const;

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
	const double targetDensity ) const;

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
	const double targetDensity ) const;

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
	const std::vector<double>& targetDensity ) const;

    /*!
     * \brief Returns a string that describes the templated
     *     EnergyPolicy.  Currently this will return either "mg" or
     *     "gray."
     */ 
    const std::string& getEnergyPolicyDescriptor() const {
	return energyPolicyDescriptor; };
};


//===========================================================================//
//                         MORE TYPEDEFS                      
//===========================================================================//

/*!
 * \brief An alias for GandolfOpacity<Gray>
 */
typedef DummyOpacity<Gray> GrayOpacity;

/*!
 * \brief An alias for GandolfOpacity<Multigroup>
 */
typedef DummyOpacity<Multigroup> MGOpacity;

/*!
 * \brief An alias for GandolfOpacity<Multigroup>
 */
typedef DummyOpacity<Multigroup> MultigroupOpacity;

} // end namespace rtt_dummy_opacity

#endif // __cdi_DummyOpacity_hh__

//---------------------------------------------------------------------------//
//                end of cdi/test/DummyOpacity.hh
//---------------------------------------------------------------------------//
