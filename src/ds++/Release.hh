//----------------------------------*-C++-*----------------------------------//
/*! 
 *  \file ds++/Release.hh
 *  \author Randy M. Roberts
 *  \date Thu Jul 15 09:31:44 1999
 *  \brief  Header file for ds++ library release function.
 */
//===========================================================================//
// $Id$
//===========================================================================//

#ifndef __ds_Release_hh__
#define __ds_Release_hh__

//===========================================================================//

#include <string>
/*!
 *
 * \brief RTT "Data Structures in C++" (DS++)  namespace.
 *
 */
namespace rtt_dsxx 
{
/*!
 * \brief  Gets the release number for the ds++ package. 
 * \return release number as a string in the form "ds-\#_\#_\#".
 * 
 * This can be used to get exact version information in codes that 
 * use ds++.
 *
 */
    const std::string release();
}

#endif                          // __ds_Release_hh__

/*!
 * \page ds++-overview Overview of the Data Structures in C++ (DS++) Package
 *
 * \version 1_0_0
 *
 * This package provides data structures and other miscellaneous support
 * for the Draco system. Classes for matrices, design-by-contract, 
 * smart-pointers, and a variety of other utilities are included.
 *
 */

//---------------------------------------------------------------------------//
//                              end of ds++/Release.hh
//---------------------------------------------------------------------------//
