//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Release.hh
 * \author Thomas M. Evans
 * \date   Thu Jul 15 09:41:09 1999
 * \brief  Release function for c4 library.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __c4_Release_hh__
#define __c4_Release_hh__

//===========================================================================//
/*!
 *
 * \page c4_overview Overview of the C4 package
 * \version 2_0_0
 * \author Tom Evans, Todd Urbatsch, Geoff Furnish

 * The c4 package is a set of classes and functions that provide wrappers for
 * message passing functionality.  Currently supported are wrappers for MPI.
 * Wrappers for SHMEM are also included; however, they are no longer fully
 * supported.
 *
 * Plans are in place to wrap the Unified Parallel Services (UPS) library.
 * The c4 package will be undergoing significant modification for release
 * 2_0_0, due out sometime in 2001.
 *
 */
//===========================================================================//
/*!
 * \namespace rtt_c4
 *
 * \brief Namespace that contains c4 package classes and variables.
 *
 * This namespace is not fully in use yet.  The only c4 services provided in
 * this namespace are the rtt_c4::release functions.  Users should use the C4
 * namespace until further notice.
 *
 */
//===========================================================================//

#include <string>

namespace rtt_c4 
{

const std::string release();

}

#endif                          // __c4_Release_hh__

//---------------------------------------------------------------------------//
//                              end of c4/Release.hh
//---------------------------------------------------------------------------//
