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

#ifndef rtt_c4_Release_hh
#define rtt_c4_Release_hh

//===========================================================================//
/*!
 *
 * \mainpage Overview of the C4 package
 * \version 2_3_0
 * \author Tom Evans, Todd Urbatsch, Geoff Furnish

 * The c4 package is a set of classes and functions that provide wrappers for
 * message passing functionality. 
 *
 * The basic functions provided by C4 sit in the C4_Functions.hh header.
 * These are canonical forms of message passing functions.  c4 is set to
 * either serial or mpi (as determined by the C4_MPI or C4_SCALAR macros).
 * To use c4 include the following
 * \code
 * #include "c4/global.hh"
 * \endcode
 */
//===========================================================================//
/*!
 * \namespace rtt_c4
 *
 * \brief Namespace that contains c4 package classes and variables.
 *
 * The C4 namespace is now deprecated and only supports pre-c4-2_0_0
 * functionality. 
 *
 */
//===========================================================================//

#include <string>

namespace rtt_c4 
{

const std::string release();

}

#endif                          // rtt_c4_Release_hh

//---------------------------------------------------------------------------//
//                              end of c4/Release.hh
//---------------------------------------------------------------------------//
