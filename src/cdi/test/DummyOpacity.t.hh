//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/GandolfOpacity.cc
 * \author Kelly Thompson
 * \date   Wed Jul 12 16:11:55 2000
 * \brief  DummyOpacity templated class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DummyOpacity.hh"

#include <cmath> // we need to define pow(double,int)

//#include "ds++/Assert.hh" // we make use of Require()

namespace rtt_dummy_opacity
{

    // ------------ //
    // Constructors //
    // ------------ //

/*!
 * \brief Constructor for DummyOpacity object.
 * 
 * See DummyOpacity.hh for details.
 *
 * Note that everything in this file must be templated by the
 * EnergyPolicy.  All Templated forms of DummyOpacity<EnergyPolicy>
 * must be instantiated in DummyOpacity_pt.cc
 */
template < class EnergyPolicy >
DummyOpacity< EnergyPolicy >::DummyOpacity( )
    : dataFilename( "none" ),
      dataDescriptor( "dummyOpacity" )
    {
	// The Multigroup class for dummyOpacity assumes that there is 
	// only one group.  If this number is changed then you must
	// manually change the multigroup class definition accordingly.
	groupBoundaries.resize(2);  
	for ( int i=0; i<groupBoundaries.size(); ++i)
	    groupBoundaries[i] = i*500.0;
    } 
 
    /*!
     * \brief Default DummyOpacity() destructor.
     *
     * This is required to correctly release memory when a
     * DummyOpacity<EnergyPolicy> is destroyed.
     */
template < class EnergyPolicy >
DummyOpacity< EnergyPolicy >::~DummyOpacity()
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
template < class EnergyPolicy >
const std::string& DummyOpacity< EnergyPolicy >::getDataDescriptor() const 
    {
	return dataDescriptor;
    }

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
     *     been done for a few STL containers in DummyOpacity_pt.cc.
     */
// template < class EnergyPolicy >
// template < class InputIterator, class OutputIterator >
// OutputIterator 
// DummyOpacity< EnergyPolicy >::getOpacity(
//     InputIterator temperatureIterator, 
//     InputIterator temperatureIteratorEnd,
//     InputIterator densityIterator, 
//     InputIterator densityIteratorEnd,
//     OutputIterator opacityIterator ) const
//     { 
// 	// call the correct accessor from the PolicyTemplate
// 	return EnergyPolicy::getOpacity(
// 	    temperatureIterator, temperatureIteratorEnd,
// 	    densityIterator, densityIteratorEnd,
// 	    opacityIterator );
//     }

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
     *     been done for a few STL containers in DummyOpacity_pt.cc.
     */
// template < class EnergyPolicy >
// template < class InputIterator, class OutputIterator >
// OutputIterator 
// DummyOpacity< EnergyPolicy >::getOpacity(
//     InputIterator temperatureIterator, 
//     InputIterator temperatureIteratorEnd,
//     const double targetDensity,
//     OutputIterator opacityIterator ) const
//     { 
// 	// call the correct accessor from the PolicyTemplate
// 	return EnergyPolicy::getOpacity(
// 	    temperatureIterator, temperatureIteratorEnd,
// 	    targetDensity,
// 	    opacityIterator );
//     }

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
     *     been done for a few STL containers in DummyOpacity_pt.cc.
     */
// template < class EnergyPolicy >
// template < class InputIterator, class OutputIterator >
// OutputIterator 
// DummyOpacity< EnergyPolicy >::getOpacity(
//     const double targetTemperature,
//     InputIterator densityIterator, 
//     InputIterator densityIteratorEnd,
//     OutputIterator opacityIterator ) const
//     { 
// 	// call the correct accessor from the PolicyTemplate
// 	return EnergyPolicy::getOpacity(
// 	    targetTemperature,
// 	    densityIterator, densityIteratorEnd,
// 	    opacityIterator );
//     }

    /*!
     * \brief Opacity accessor that returns a single opacity (or a
     *     vector of opacities for the multigroup EnergyPolicy) that 
     *     corresponds to the provided temperature and density.
     */
template < class EnergyPolicy >
DummyOpacity< EnergyPolicy >::OpacityType 
DummyOpacity< EnergyPolicy >::getOpacity(
    const double targetTemperature,
    const double targetDensity ) const
    { 
	// call the correct accessor from the PolicyTemplate
	return EnergyPolicy::getOpacity( targetTemperature,
					 targetDensity ); 
    }

    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     temperatures and a single density value.
     */
template < class EnergyPolicy >
std::vector< typename DummyOpacity< EnergyPolicy >::OpacityType > 
DummyOpacity< EnergyPolicy >::getOpacity(
    const std::vector<double>& targetTemperature,
    const double targetDensity ) const
    { 
	// call the correct accessor from the PolicyTemplate
	return EnergyPolicy::getOpacity( targetTemperature,
					 targetDensity );
    }

    /*!
     * \brief Opacity accessor that returns a vector of opacities (or a
     *     vector of vectors of opacities for the multigroup
     *     EnergyPolicy) that correspond to the provided vector of
     *     densities and a single temperature value.
     */
template < class EnergyPolicy >
std::vector< typename DummyOpacity< EnergyPolicy >::OpacityType > 
DummyOpacity< EnergyPolicy >::getOpacity(
    const double targetTemperature,
    const std::vector<double>& targetDensity ) const
    { 
	// call the correct accessor from the PolicyTemplate
	return EnergyPolicy::getOpacity( targetTemperature,
					 targetDensity );
    }

    // It is not clear how to assume order of opacity(temp,dens) when
    // accessed in this manner --> for now use the STL-style accessor
    // or a loop over one of the other vector-accessors.

// template < class EnergyPolicy >
// std::vector< typename DummyOpacity< EnergyPolicy >::OpacityType > 
// DummyOpacity< EnergyPolicy >::getOpacity(
//     const std::vector<double>& targetTemperature,
//     const std::vector<double>& targetDensity ) const
//     { 
// 	return EnergyPolicy::getOpacity( targetTemperature,
// 					 targetDensity,
// 					 spGandolfDataTable ); 
//     }

    /*!
     * \brief Returns a vector of temperatures that define the cached
     *     opacity data table.
     */
template < class EnergyPolicy >
std::vector<double> DummyOpacity< EnergyPolicy >::getTemperatureGrid() const
    {
	// Retrieve the temperature grid for the Data Object.  These
	// are log(temp) values.
	std::vector<double> tGrid(2);
	// convert the log(temp) grid to a temp grid.
	for ( int i=0; i<tGrid.size(); ++i )
	    tGrid[i] = (i+1)*0.1;
	// return the temperature grid.
	return tGrid;
    }

    /*!
     * \brief Returns the size of the temperature grid.
     */
template < class EnergyPolicy >
int DummyOpacity< EnergyPolicy >::getNumTemperatures() const
    {
	return 2;
    }

    /*!
     * \brief Returns a vector of densities that define the cached
     *     opacity data table.
     */
template < class EnergyPolicy >
std::vector<double> DummyOpacity< EnergyPolicy >::getDensityGrid() const
    {
	// Retrieve the density grid for the Data Object.  These
	// are log(rho) values.
	std::vector<double> rhoGrid(3);
	// convert the log(rho) grid to a density grid.
	for ( int i=0; i<rhoGrid.size(); ++i )
	    rhoGrid[i] = (i+1)*0.005;
	// return the density grid.
	return rhoGrid;
    }

    /*! 
     * \brief Returns the size of the density grid.
     */
template < class EnergyPolicy >
int DummyOpacity< EnergyPolicy >::getNumDensities() const
    {
	return 3;
    }

    /*!
     * \brief Returns a vector of energy values (keV) that define the
     *     energy boundaries of the cached multigroup opacity data
     *     table.  (This accessor is only valid for the Multigroup
     *     EnergyPolicy version of DummyOpacity.)
     */
template < class EnergyPolicy >
const std::vector<double>& DummyOpacity< EnergyPolicy >::getGroupBoundaries() const
    {
	return groupBoundaries;
    }

    /*!
     * \brief Returns the number of group boundaries found in the
     *     current multigroup data set.
     */
template < class EnergyPolicy >
int DummyOpacity< EnergyPolicy >::getNumGroupBoundaries() const
    {
	return groupBoundaries.size();
    }

		       // ----------------- //
		       // Gray Policy Class //
		       // ----------------- //

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
     *     been done for a few STL containers in DummyOpacity_pt.cc.
     */
// template < class InputIterator, class OutputIterator >
// OutputIterator Gray::getOpacity(
//     InputIterator temperatureIterator, 
//     InputIterator temperatureIteratorEnd,
//     InputIterator densityIterator, 
//     InputIterator densityIteratorEnd,
//     OutputIterator opacityIterator ) const
//     {
// 	// Loop over all (temperature,density) tuple values.
// 	for ( ; temperatureIterator != temperatureIteratorEnd;
// 	      ++temperatureIterator, ++densityIterator, ++opacityIterator )
// 	    *opacityIterator = *densityIterator
// 		* pow( *temperatureIterator, 4 );
// 	return opacityIterator;
//     }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of temperatures in an STL container.
     *     An set of opacity values will be returned for each temperature
     *     provided.  The opacity iterator should have length equal to 
     *     the length of the temperatureIterator.
     */
// template < class InputIterator, class OutputIterator >
// OutputIterator Gray::getOpacity(
//     InputIterator temperatureIterator, 
//     InputIterator temperatureIteratorEnd,
//     const double targetDensity,
//     OutputIterator opacityIterator ) const
//     {
// 	for ( ; temperatureIterator != temperatureIteratorEnd;
// 	      ++temperatureIterator, ++opacityIterator )
// 	    *opacityIterator = targetDensity
// 		* pow( *temperatureIterator, 4 );
// 	return opacityIterator;
//     }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of densities in an STL container.
     *     An opacity value will be returned for each density
     *     provided.  The opacity iterator should have length equal to 
     *     to the length of the densityIterator.
     */
// template < class InputIterator, class OutputIterator >
// OutputIterator Gray::getOpacity(
//     const double targetTemperature,
//     InputIterator densityIterator, 
//     InputIterator densityIteratorEnd,
//     OutputIterator opacityIterator ) const
//     {
// 	for ( ; densityIterator != densityIteratorEnd;
// 	      ++densityIterator, ++opacityIterator )
// 	    *opacityIterator = *densityIterator
// 		* pow( targetTemperature, 4 );
// 	return opacityIterator;
//     }

   /*!
     * \brief Opacity accessor that returns a single opacity that 
     *     corresponds to the provided temperature and density.
     */
Gray::OpacityType Gray::getOpacity( 
    const double targetTemperature,
    const double targetDensity ) const
    {
	return targetDensity * pow( targetTemperature, 4 );
    }

    /*!
     * \brief Opacity accessor that returns a vector of
     *     opacities that correspond to the provided vector of
     *     temperatures and a single density value.
     */
std::vector< Gray::OpacityType > Gray::getOpacity( 
    const std::vector<double>& targetTemperature,
    const double targetDensity ) const
    {
	std::vector< OpacityType > opacity( targetTemperature.size() );
	for ( int i=0; i<targetTemperature.size(); ++i )
	    opacity[i] = targetDensity * pow( targetTemperature[i], 4 );
	return opacity;
    }

    /*!
     * \brief Opacity accessor that returns a vector of
     *     opacities that correspond to the provided vector of
     *     densities and a single temperature value.
     */
std::vector< Gray::OpacityType > Gray::getOpacity( 
    const double targetTemperature,
    const std::vector<double>& targetDensity ) const
    {
	std::vector< OpacityType > opacity( targetDensity.size() );
	for ( int i=0; i<targetDensity.size(); ++i )
	    opacity[i] = targetDensity[i] * pow( targetTemperature, 4 );
	return opacity;
    }


		    // ----------------------- //
		    // Multigroup Policy Class //
		    // ----------------------- //

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
     *     been done for a few STL containers in DummyOpacity_pt.cc.
     */
// template < class InputIterator, class OutputIterator >
// OutputIterator Multigroup::getOpacity(
//     InputIterator temperatureIterator, 
//     InputIterator temperatureIteratorEnd,
//     InputIterator densityIterator, 
//     InputIterator densityIteratorEnd,
//     OutputIterator opacityIterator ) const
//     {
// 	// number of groups in this multigroup set.
// 	const int ng = 1;  
	
// 	// loop over the (temperature,density) tuple.
// 	for ( ; temperatureIterator != temperatureIteratorEnd;
// 	      ++temperatureIterator, ++densityIterator )
// 	    for ( int i=0; i<ng; ++i, ++opacityIterator )
// 		*opacityIterator = *densityIterator
// 		    * pow( *temperatureIterator, 4 );

// 	return opacityIterator;
//     }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of temperatures in an STL container.
     *     An set of opacity values will be returned for each temperature
     *     provided.  The opacity iterator should have length equal to 
     *     the length of the temperatureIterator times the number of
     *     energy groups.
     */
// template < class InputIterator, class OutputIterator >
// OutputIterator Multigroup::getOpacity(
//     const double targetTemperature,
//     InputIterator densityIterator, 
//     InputIterator densityIteratorEnd,
//     OutputIterator opacityIterator ) const
//     {
// 	const int ng = 1;
// 	for ( ; densityIterator != densityIteratorEnd; ++densityIterator )
// 	    for ( int i=0; i<ng; ++i, ++opacityIterator )
// 		*opacityIterator = *densityIterator
// 		    * pow( targetTemperature, 4 );

// 	return opacityIterator;
//     }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of densities in an STL container.
     *     An opacity value will be returned for each density
     *     provided.  The opacity iterator should have length equal to 
     *     to the length of the densityIterator times the number of
     *     energy groups.
     */
// template < class InputIterator, class OutputIterator >
// OutputIterator Multigroup::getOpacity(
//     InputIterator temperatureIterator, 
//     InputIterator temperatureIteratorEnd,
//     const double targetDensity,
//     OutputIterator opacityIterator ) const
//     {
// 	const int ng = 1;
// 	for ( ; temperatureIterator != temperatureIteratorEnd;
// 	     ++temperatureIterator )
// 	    for ( int i=0; i<ng; ++i, ++opacityIterator )
// 		*opacityIterator = targetDensity
// 		    * pow( *temperatureIterator, 4 );		

// 	return opacityIterator;
//     }

    /*!
     * \brief Opacity accessor that returns a vector of opacities that 
     *     corresponds to the provided temperature and density.
     */
Multigroup::OpacityType Multigroup::getOpacity(
    const double targetTemperature,
    const double targetDensity ) const
    {
	int numGroups = 1;
	OpacityType opacity( numGroups );
	for ( int i=0; i<numGroups; ++i)
	    opacity[i] = targetDensity * pow( targetTemperature, 4 );
	return opacity;
    }

    /*!
     * \brief Opacity accessor that returns a vector of vectors of
     *     opacities that correspond to the provided vector of
     *     temperatures and a single density value.
     */
std::vector< Multigroup::OpacityType > Multigroup::getOpacity(
    const std::vector<double>& targetTemperature,
    const double targetDensity ) const
    {
	int ng = 1;
	std::vector< OpacityType > opacity( targetTemperature.size() );
	for ( int i=0; i<targetTemperature.size(); ++i )
	    {
		opacity[i].resize(ng);
		for ( int ig=0; ig<ng; ++ig)
		    opacity[i][ig] = targetDensity * pow( targetTemperature[i], 4 );
	    }
	return opacity;
    }

    /*!
     * \brief Opacity accessor that returns a vector of vectors of
     *     opacities that correspond to the provided vector of
     *     densities and a single temperature value.
     */
std::vector< Multigroup::OpacityType > Multigroup::getOpacity(
    const double targetTemperature,
    const std::vector<double>& targetDensity ) const
    {
	int ng = 1;
	std::vector< OpacityType > opacity( targetDensity.size() );
	for ( int i=0; i<targetDensity.size(); ++i )
	    {
		opacity[i].resize(ng);
		for ( int ig=0; ig<ng; ++ig)
		    opacity[i][ig] = targetDensity[i] * pow( targetTemperature, 4 );
	    }
	return opacity;
    }


} // end namespace rtt_dummy_opacity


//---------------------------------------------------------------------------//
//                              end of DummyOpacity.cc
//---------------------------------------------------------------------------//
