//----------------------------------*-C++-*----------------------------------//
/*! 
 *  \file   ds++/Release.hh
 *  \author Thomas Evans
 *  \date   Thu Jul 15 09:31:44 1999
 *  \brief  Header file for ds++ library release function.
 */
//===========================================================================//
// $Id$
//===========================================================================//

#ifndef rtt_ds_Release_hh
#define rtt_ds_Release_hh

#include <string>

//===========================================================================//
/*!
 * \mainpage Overview of the Data Structures in C++ (DS++) Package
 * \author Tom Evans, Kelly Thompson, Todd Urbatsch
 * \version 1_9_0
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

#endif                          // rtt_ds_Release_hh

//---------------------------------------------------------------------------//
//                              end of ds++/Release.hh
//---------------------------------------------------------------------------//
