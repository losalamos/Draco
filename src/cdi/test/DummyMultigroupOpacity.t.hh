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
     * \brief Opacity accessor that returns an STL container of
     *     opacities that correspond to a tuple of provided STL
     *     containers (temperatures and densities).  The length of the 
     *     temperature and the the density container should be equal
     *     and the length of the opacity container should be
     *     numGroups x temperature.size().
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
     * \brief Opacity accessor that returns an STL container of
     *     opacities that correspond to a list of provided STL
     *     temperature values.  The length of the opacity container
     *     should be numGroups x temperature.size().
     */
    template < class TemperatureIterator, class OpacityIterator >
	OpacityIterator DummyMultigroupOpacity::getOpacity(
	    TemperatureIterator tempIter,
    TemperatureIterator templast,
	    double targetDensity,
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
     * \brief Opacity accessor that returns an STL container of
     *     opacities that correspond to a list of provided STL
     *     density values and a fixed temperature.  The length of the
     *     opacity container should be numGroups x density.size().
     */
    template < class DensityIterator, class OpacityIterator >
	OpacityIterator DummyMultigroupOpacity::getOpacity(
	    double targetTemperature,
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
