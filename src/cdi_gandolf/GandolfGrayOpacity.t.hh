//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfGrayOpacity.t.hh
 * \author Kelly Thompson
 * \date   Wed Jan 24 9:44:55 2001
 * \brief  GandolfGrayOpacity templated class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfGrayOpacity_t_hh__
#define __cdi_gandolf_GandolfGrayOpacity_t_hh__

#include "GandolfGrayOpacity.hh"

#include "GandolfWrapper.hh"    // we make calls to the wrapper routines.

#include "GandolfDataTable.hh"  // we have a smart pointer to a
                                // GandolfDataTable object.

#include <cmath> // we need to define log(double) and exp(double)

#include "ds++/Assert.hh" // we make use of Require()

namespace rtt_cdi_gandolf
{

    // --------------------------------- //
    // STL-like accessors for getOpacity //
    // --------------------------------- //

    // ------------------------------------------ //
    // getOpacity with Tuple of (T,rho) arguments //
    // ------------------------------------------ //

template < class TemperatureIterator, class DensityIterator,
           class OpacityIterator >
OpacityIterator GandolfGrayOpacity::getOpacity(
    TemperatureIterator tempIter, 
    TemperatureIterator tempLast,
    DensityIterator densIter, 
    DensityIterator densLast,
    OpacityIterator opIter ) const
    { 
	// from twix:/scratch/tme/kai/KCC_BASE/include/algorithm

	// assert that the two input iterators have compatible sizes.
 	Require( std::distance( tempIter, tempLast )
 		 == std::distance( densIter, densLast ) );

	// Loop over all (temperature,density) tuple values.
	for ( ; tempIter != tempLast;
	      ++tempIter, ++densIter, ++opIter )
	    // Call the Gandolf Logorithmic Interpolator for Gray data.
	    wrapper::wgintgrlog( spGandolfDataTable->getLogTemperatures(),
				 spGandolfDataTable->getNumTemperatures(), 
				 spGandolfDataTable->getLogDensities(), 
				 spGandolfDataTable->getNumDensities(),
				 spGandolfDataTable->getLogOpacities(), 
				 spGandolfDataTable->getNumOpacities(),
				 log( *tempIter ),
				 log( *densIter ), 
				 *opIter );
	return opIter;
    }

// ------------------------------------ // 
// getOpacity() with container of temps //
// ------------------------------------ // 

template < class TemperatureIterator, class OpacityIterator >
OpacityIterator GandolfGrayOpacity::getOpacity(
    TemperatureIterator tempIter,
    TemperatureIterator tempLast,
    double targetDensity,
    OpacityIterator opIter ) const
    { 
	// loop over all the entries the temperature container and
	// calculate an opacity value for each.
	for ( ; tempIter != tempLast; ++tempIter, ++opIter )
	    // Call the Gandolf Logorithmic Interpolator for Gray data.
	    wrapper::wgintgrlog( spGandolfDataTable->getLogTemperatures(),
				 spGandolfDataTable->getNumTemperatures(), 
				 spGandolfDataTable->getLogDensities(), 
				 spGandolfDataTable->getNumDensities(),
				 spGandolfDataTable->getLogOpacities(), 
				 spGandolfDataTable->getNumOpacities(),
				 log( *tempIter ),
				 log( targetDensity ), 
				 *opIter );
	return opIter;
    }

// ---------------------------------------- // 
// getOpacity() with container of densities //
// ---------------------------------------- //

template < class DensityIterator, class OpacityIterator >
OpacityIterator GandolfGrayOpacity::getOpacity(
    double targetTemperature,
    DensityIterator densIter, 
    DensityIterator densLast,
    OpacityIterator opIter ) const
    { 
	// loop over all the entries the density container and
	// calculate an opacity value for each.
	for ( ; densIter != densLast; ++densIter, ++opIter )
	    // Call the Gandolf Logorithmic Interpolator for Gray data.
	    wrapper::wgintgrlog( spGandolfDataTable->getLogTemperatures(),
				 spGandolfDataTable->getNumTemperatures(), 
				 spGandolfDataTable->getLogDensities(), 
				 spGandolfDataTable->getNumDensities(),
				 spGandolfDataTable->getLogOpacities(), 
				 spGandolfDataTable->getNumOpacities(),
				 log( targetTemperature ),
				 log( *densIter ), 
				 *opIter );
	return opIter;
    }

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfGrayOpacity_t_hh__

//---------------------------------------------------------------------------//
// end of GandolfGrayOpacity.t.hh
//---------------------------------------------------------------------------//
