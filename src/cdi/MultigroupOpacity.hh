//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/MultigroupOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 8 14:58:55 2001
 * \brief  MultigroupOpacity class header file (an abstract class)
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_MultigroupOpacity_hh__
#define __cdi_MultigroupOpacity_hh__

#include <vector>
//#include <string>

//#include "ds++/SP.hh"

namespace rtt_cdiMultigroupOpacity
{
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
 * \class MultigroupOpacity
 *
 * \brief 
 *
 * \sa  
 */

/*!
 * \example cdi/test/tDummyOpacity.cc
 * \example cdi/test/tCDI.cc
 *
 */
//===========================================================================//

class MultigroupOpacity
{
    // DATA

  public:

    // ----------- //
    // Destructors //
    // ----------- //

    /*!
     * \brief Default Opacity() destructor.
     *
     * This is required to correctly release memory when an
     * Opacity<EnergyPolicy> is destroyed.
     */
     virtual ~MultigroupOpacity() {};

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
//     template < class InputIterator, class OutputIterator >
//     virtual OutputIterator getOpacity( 
// 	InputIterator temperatureIterator, 
// 	InputIterator temperatureIteratorEnd,
// 	InputIterator densityIterator, 
// 	InputIterator densityIteratorEnd,
// 	OutputIterator opacityIterator ) const = 0;

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
//     template < class InputIterator, class OutputIterator >
//     virtual OutputIterator getOpacity( 
// 	InputIterator temperatureIterator,
// 	InputIterator temperatureIteratorEnd,
// 	const double targetDensity,
// 	OutputIterator opacityIterator ) const = 0;

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
//     template < class InputIterator, class OutputIterator >
//     virtual OutputIterator getOpacity( 
// 	const double targetTemperature,
// 	InputIterator densityIterator, 
// 	InputIterator densityIteratorEnd,
// 	OutputIterator opacityIterator ) const = 0;

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
    virtual std::vector<double> getOpacity( 
	const double targetTemperature,
	const double targetDensity ) const = 0; 

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
    virtual std::vector< std::vector<double> > getOpacity( 
	const std::vector<double>& targetTemperature,
	const double targetDensity ) const = 0; 
    
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
    virtual std::vector< std::vector<double> > getOpacity( 
	const double targetTemperature,
	const std::vector<double>& targetDensity ) const = 0; 

    /*!
     * \brief Returns a string that describes the templated
     *     EnergyPolicy.  Currently this will return either "mg" or
     *     "gray."
     */ 
    virtual const std::string& getEnergyPolicyDescriptor() const = 0;

    /*!
     * \brief Returns a "plain English" description of the opacity
     *     data that this class references. (e.g. "Gray Rosseland
     *     Scattering".) 
     *
     * The definition of this function is not included here to prevent 
     *     the inclusion of the GandolfFile.hh definitions within this 
     *     header file.
     */
    virtual const std::string& getDataDescriptor() const = 0;

    /*!
     * \brief Returns the name of the associated IPCRESS file.
     *
     * The definition of this function is not included here to prevent 
     *     the inclusion of the GandolfFile.hh definitions within this 
     *     header file.
     */
    virtual const std::string& getDataFilename() const = 0;

    /*!
     * \brief Returns a vector of temperatures that define the cached
     *     opacity data table.
     * 
     * We do not return a const reference because this function
     * must construct this information from more fundamental tables.
     */
    virtual std::vector<double> getTemperatureGrid() const = 0;

    /*!
     * \brief Returns a vector of densities that define the cached
     *     opacity data table.
     * 
     * We do not return a const reference because this function
     * must construct this information from more fundamental tables.
     */
    virtual std::vector<double> getDensityGrid() const = 0;

    /*!
     * \brief Returns a vector of energy values (keV) that define the
     *     energy boundaries of the cached multigroup opacity data
     *     table.  (This accessor is only valid for the Multigroup
     *     EnergyPolicy version of GandolfOpacity.)
     */
    virtual const std::vector<double>& getGroupBoundaries() const = 0;
    
    /*!
     * \brief Returns the size of the temperature grid.
     */
    virtual int getNumTemperatures() const = 0;

    /*! 
     * \brief Returns the size of the density grid.
     */
    virtual int getNumDensities() const = 0;

    /*!
     * \brief Returns the number of group boundaries found in the
     *     current multigroup data set.
     */
    virtual int getNumGroupBoundaries() const = 0;

}; // end of class MultigroupOpacity


} // end namespace rtt_cdi_multigroup_opacity

#endif // __cdi_MultigroupOpacity_hh__

//---------------------------------------------------------------------------//
//                end of cdi/MultigroupOpacity.hh
//---------------------------------------------------------------------------//
