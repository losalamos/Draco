//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   pcgWrap/Release.hh
 * \author Rob Lowrie and Tom Evans
 * \date   Thu Jul 15 09:44:12 1999
 * \brief  Release function for pcgWrap library.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __pcgWrap_Release_hh__
#define __pcgWrap_Release_hh__

//===========================================================================//
/*!
 * \page pcgWrap_overview Overview of the pcgWrap package
 * \version 1_2_0
 * \author Rob Lowrie and Dave Nystrom
 *
 * pcgWrap is a rather messy interface to the Preconditioned Conjugate
 * Gradient library (pcg).  Because the pcg package is no longer actively
 * maintained, pcgWrap is an essentially "dead" package.  It will be
 * superceded by interfaces to the UPS linear solver packages.
 *
 * Needless to say, to use pcgWrap, libpcg must be linked.  This can be added
 * for testing purposes during configure with the --with-pcg and
 * --with-pcg-lib options.  Also, pcg lib depends on lapack, for non-default
 * behavior use the --with-lapack and --with-lapack-lib options to
 * configure. 
 *
 * Because this package does not have long for this world it is not
 * extensively documented.  For guidance in its use (if you need it
 * temporarily) see R. Lowrie <lowrie@lanl.gov>.
 *
 */
//===========================================================================//
/*!
 * \namespace rtt_pcgWrap
 *
 * \brief Namespace that contains the pcgWrap package classes and variables.
 */
//===========================================================================//

#include <string>

namespace rtt_pcgWrap
{
    const std::string release();
}

#endif                          // __pcgWrap_Release_hh__

//---------------------------------------------------------------------------//
//                              end of pcgWrap/Release.hh
//---------------------------------------------------------------------------//
