//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyGrayOpacity.t.hh
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:33:51 2001
 * \brief  DummyGrayOpacity class templated implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_DummyGrayOpacity_t_hh__
#define __cdi_DummyGrayOpacity_t_hh__

#include "DummyGrayOpacity.hh"
#include <vector>

namespace rtt_dummyGrayOpacity
{

    // --------- //
    // Accessors //
    // --------- //
    
template < class OpacityIterator, class TemperatureIterator >
OpacityIterator DummyGrayOpacity::getOpacity( 
    TemperatureIterator tempIter,
    TemperatureIterator tempLast,
    const double targetDensity,
    OpacityIterator opacityIter ) const
    {
	// loop over all temperatures in the range
	// (tempFirst,tempLast).
	for ( ; tempIter != tempLast; ++tempIter, ++opacityIter )
	    *opacityIter = *tempIter + targetDensity/1000.0;
	return opacityIter;
    }

template < class OpacityIterator, class DensityIterator >
OpacityIterator DummyGrayOpacity::getOpacity( 
    const double targetTemperature,
    DensityIterator densIter,
    DensityIterator densLast,
    OpacityIterator opacityIter ) const
    {
	// loop over all temperatures in the range
	// (tempFirst,tempLast).
	for ( ; densIter != densLast; ++densIter, ++opacityIter )
	    *opacityIter = targetTemperature + *densIter/1000.0;
	return opacityIter;
    }

template < class OpacityIterator, class TemperatureIterator,
           class DensityIterator >
OpacityIterator DummyGrayOpacity::getOpacity( 
    TemperatureIterator tempIter,
    TemperatureIterator tempLast,
    DensityIterator densIter,
    DensityIterator densLast,
    OpacityIterator opacityIter ) const
    {
	// loop over all temperatures in the range
	// (tempFirst,tempLast).
	for ( ; densIter != densLast; ++tempIter, ++densIter, ++opacityIter )
	    *opacityIter = *tempIter + *densIter/1000.0;
	return opacityIter;
    }

} // end rtt_dummyGrayOpacity

#endif // __cdi_test_DummyGrayOpacity_t_hh__

//---------------------------------------------------------------------------//
// end of DummyGrayOpacity.cc
//---------------------------------------------------------------------------//
