//----------------------------------*-C++-*----------------------------------//
/*! 
 *  \file   ds++/Release.hh
 *  \author Randy M. Roberts
 *  \date   Thu Jul 15 09:31:44 1999
 *  \brief  Header file for ds++ library release function.
 */
//===========================================================================//
// $Id$
//===========================================================================//

#ifndef __ds_Release_hh__
#define __ds_Release_hh__

#include <string>

//===========================================================================//
/*!
 * \page ds++-overview Overview of the Data Structures in C++ (DS++) Package
 *
 * \version 1_4_0
 *
 * This package provides data structures and other miscellaneous support
 * for the Draco system. Classes for matrices, design-by-contract, 
 * smart-pointers, and a variety of other utilities are included.
 *
 */
//===========================================================================//
/*!
 * \namespace rtt_dsxx
 *
 * \brief RTT "Data Structures in C++" (DS++)  namespace.
 *
 */
//===========================================================================//

namespace rtt_dsxx 
{

//! Query package for the release number.
const std::string release();

} // end of rtt_ds++

#endif                          // __ds_Release_hh__

//---------------------------------------------------------------------------//
//                              end of ds++/Release.hh
//---------------------------------------------------------------------------//
