//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   meshReaders/Release.hh
 * \author B.T. Adams
 * \date   Fri Aug 27 10:33:26 1999
 * \brief  Header file for meshReaders library release function.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshReaders_Release_hh__
#define __meshReaders_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of meshReaders; 
// this can be used to get exact version information in codes that 
// use meshReaders
// 
//===========================================================================//

#include <string>

namespace rtt_format 
{
/*!
 * \brief  Gets the release number for the meshReaders package. 
 * \return release number as a string in the form "meshReaders-\#_\#_\#"
 */
    const std::string release();
}

#endif                          // __meshReaders_Release_hh__



//---------------------------------------------------------------------------//
//                              end of Release.hh
//---------------------------------------------------------------------------//
