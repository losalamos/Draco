//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/GandolfDummyOpacity.t.hh
 * \author Kelly Thompson
 * \date   Wed Jan 17 14:34:16 2001
 * \brief  DummyMultigroupOpacity class templeted implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_DummyMultigroupOpacity_t_hh__
#define __cdi_DummyMultigroupOpacity_t_hh__

#include "DummyMultigroupOpacity.hh"
#include <vector>

namespace rtt_dummyMultigroupOpacity
{

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
     *     been done for a few STL containers in DummyMultigroupOpacity_pt.cc.
     */
template < class TemperatureIterator, class DensityIterator, 
           class OpacityIterator >
OpacityIterator DummyMultigroupOpacity::getOpacity(
    TemperatureIterator tempIter,
    TemperatureIterator tempLast,
    DensityIterator densIter,
    DensityIterator densLast,
    OpacityIterator opacityIter ) const
    { 
	int ng = numGroupBoundaries - 1;
	// loop over all temperatures and densities in the range
	// (tempFirst,tempLast) & (densIter,densLast).
	for ( ; densIter != densLast && tempIter != tempLast; ++tempIter, ++densIter )
	    for ( int ig=0; ig<ng; ++ig, ++opacityIter )
		*opacityIter = 2.0 * ( *tempIter + *densIter/1000.0 )
 		    / ( groupBoundaries[ig] + groupBoundaries[ig+1] );
	return opacityIter;
    }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of temperatures in an STL container.
     *     An opacity value will be returned for each temperature
     *     provided.  The opacity iterator should be the same length
     *     as the temperatureIterator for gray data or the length of
     *     the temperatureIterator times the number of energy groups
     *     for multigroup data.
     *
     * The InputIterator and OpacityIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in DummyMultigroupOpacity_pt.cc.
     */
template < class TemperatureIterator, class OpacityIterator >
OpacityIterator DummyMultigroupOpacity::getOpacity(
    TemperatureIterator tempIter,
    TemperatureIterator templast,
    const double targetDensity,
    OpacityIterator opacityIter ) const
    { 
	int ng = numGroupBoundaries - 1;
	// loop over all temperatures in the range
	// (tempFirst,tempLast).
	for ( ; tempIter != tempLast; ++tempIter )
	    for ( int ig=0; ig<ng; ++ig, ++opacityIter )
		*opacityIter = 2.0 * ( *tempIter + targetDensity/1000.0 )
		    / ( groupBoundaries[ig] + groupBoundaries[ig+1] );
	return opacityIter;
    }

    /*!
     * \brief Opacity accessor that utilizes STL-like iterators.  This 
     *     accessor expects a list of densities in an STL container.
     *     An opacity value will be returned for each density
     *     provided.  The opacity iterator should be the same length
     *     as the densityIterator for gray data or the length of the
     *     densityIterator times the number of energy groups for
     *     multigroup data.
     *
     * The InputIterator and OpacityIterator classes must be
     *     instantiated for each STL container used.  This has already
     *     been done for a few STL containers in DummyMultigroupOpacity_pt.cc.
     */
template < class DensityIterator, class OpacityIterator >
OpacityIterator DummyMultigroupOpacity::getOpacity(
    const double targetTemperature,
    DensityIterator densIter,
    DensityIterator densLast,
    OpacityIterator opacityIter ) const
    { 
	int ng = numGroupBoundaries - 1;
	// loop over all densities in the range
	// (densIter,densLast).
	for ( ; densIter != densLast; ++densIter )
	    for ( int ig=0; ig<ng; ++ig, ++opacityIter )
		*opacityIter = 2.0 * ( targetTemperature + *densIter/1000.0 )
		    / ( groupBoundaries[ig] + groupBoundaries[ig+1] );
	return opacityIter;
    }

} // end namespace rtt_dummyMultigroupOpacity

#endif // __cdi_DummyMultigroupOpacity_t_hh__

//---------------------------------------------------------------------------//
// end of DummyMultigroupOpacity.t.hh
//---------------------------------------------------------------------------//
