//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfMultigroupOpacity.t.hh
 * \author Kelly Thompson
 * \date   Thu Jan 25 13:53:65 2001
 * \brief  GandolfMultigroupOpacity templated class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfMultigroupOpacity_t_hh__
#define __cdi_gandolf_GandolfMultigroupOpacity_t_hh__

#include "GandolfMultigroupOpacity.hh"

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
OpacityIterator GandolfMultigroupOpacity::getOpacity(
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

	// number of groups in this multigroup set.
	const int ng = spGandolfDataTable->getNumGroupBoundaries()-1;
	
	// temporary opacity vector used by the wrapper.  The returned 
	// data will be copied into the opacityIterator.
	std::vector<double> mgOpacity( ng );

	// loop over the (temperature,density) tuple.
	for ( ; tempIter != tempLast; ++tempIter, ++densIter )
	    {
		// Call the Gandolf interpolator.
		// the vector opacity is returned.
		wrapper::wgintmglog( spGandolfDataTable->getLogTemperatures(),
				     spGandolfDataTable->getNumTemperatures(), 
				     spGandolfDataTable->getLogDensities(), 
				     spGandolfDataTable->getNumDensities(),
				     spGandolfDataTable->getNumGroupBoundaries(),
				     spGandolfDataTable->getLogOpacities(), 
				     spGandolfDataTable->getNumOpacities(),
				     log( *tempIter ),
				     log( *densIter ), 
				     mgOpacity );
		
		// The opacity vector contains the solution.  Now
		// we copy this solution into the OpacityIterator
		for ( int i=0; i<ng; ++i, ++opIter )
		    *opIter = mgOpacity[i];
	    }
	return opIter;
    }

// ------------------------------------ // 
// getOpacity() with container of temps //
// ------------------------------------ // 

template < class TemperatureIterator, class OpacityIterator >
OpacityIterator GandolfMultigroupOpacity::getOpacity(
    TemperatureIterator tempIter,
    TemperatureIterator tempLast,
    const double targetDensity,
    OpacityIterator opIter ) const
    { 
	// number of groups in this multigroup set.
	const int ng = spGandolfDataTable->getNumGroupBoundaries()-1;
	
	// temporary opacity vector used by the wrapper.  The returned 
	// data will be copied into the opacityIterator.
	std::vector<double> mgOpacity( ng );

	// loop over the (temperature,density) tuple.
	for ( ; tempIter != tempLast; ++tempIter )
	    {
		// Call the Gandolf interpolator.
		// the vector opacity is returned.
		wrapper::wgintmglog( spGandolfDataTable->getLogTemperatures(),
				     spGandolfDataTable->getNumTemperatures(), 
				     spGandolfDataTable->getLogDensities(), 
				     spGandolfDataTable->getNumDensities(),
				     spGandolfDataTable->getNumGroupBoundaries(),
				     spGandolfDataTable->getLogOpacities(), 
				     spGandolfDataTable->getNumOpacities(),
				     log( *tempIter ),
				     log( targetDensity ), 
				     mgOpacity );
		
		// The opacity vector contains the solution.  Now
		// we copy this solution into the OpacityIterator
		for ( int i=0; i<ng; ++i, ++opIter )
		    *opIter = mgOpacity[i];
	    }
	return opIter;
    }

// ---------------------------------------- // 
// getOpacity() with container of densities //
// ---------------------------------------- //

template < class DensityIterator, class OpacityIterator >
OpacityIterator GandolfMultigroupOpacity::getOpacity(
    const double targetTemperature,
    DensityIterator densIter, 
    DensityIterator densLast,
    OpacityIterator opIter ) const
    { 
	// number of groups in this multigroup set.
	const int ng = spGandolfDataTable->getNumGroupBoundaries()-1;
	
	// temporary opacity vector used by the wrapper.  The returned 
	// data will be copied into the opacityIterator.
	std::vector<double> mgOpacity( ng );

	// loop over the (temperature,density) tuple.
	for ( ; densIter != densLast; ++densIter )
	    {
		// Call the Gandolf interpolator.
		// the vector opacity is returned.
		wrapper::wgintmglog( spGandolfDataTable->getLogTemperatures(),
				     spGandolfDataTable->getNumTemperatures(), 
				     spGandolfDataTable->getLogDensities(), 
				     spGandolfDataTable->getNumDensities(),
				     spGandolfDataTable->getNumGroupBoundaries(),
				     spGandolfDataTable->getLogOpacities(), 
				     spGandolfDataTable->getNumOpacities(),
				     log( targetTemperature ),
				     log( *densIter ), 
				     mgOpacity );
		
		// The opacity vector contains the solution.  Now
		// we copy this solution into the OpacityIterator
		for ( int i=0; i<ng; ++i, ++opIter )
		    *opIter = mgOpacity[i];
	    }
	return opIter;
    }

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfMultigroupOpacity_t_hh__

//---------------------------------------------------------------------------//
// end of GandolfMultigroupOpacity.t.hh
//---------------------------------------------------------------------------//


