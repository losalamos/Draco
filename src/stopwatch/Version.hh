//----------------------------------*-C++-*----------------------------------//
// Version.hh
// Shawn Pautz
// Tue Mar 23 15:37:40 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Version function for stopwatch library
//---------------------------------------------------------------------------//

#ifndef __stopwatch_Version_hh__
#define __stopwatch_Version_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of stopwatch; 
// this can be used to get exact version information in codes that 
// use stopwatch
// 
//===========================================================================//

#include <string>

namespace rtt_stopwatch 
{
    const std::string version();
}

#endif                          // __stopwatch_Version_hh__

//---------------------------------------------------------------------------//
//                              end of Version.hh
//---------------------------------------------------------------------------//
