//----------------------------------*-C++-*----------------------------------//
// Release.hh
// Randy M. Roberts
// Thu Jul 15 09:44:12 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function for pcgWrap library
//---------------------------------------------------------------------------//

#ifndef __pcgWrap_Release_hh__
#define __pcgWrap_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of pcgWrap; 
// this can be used to get exact version information in codes that 
// use pcgWrap
// 
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
