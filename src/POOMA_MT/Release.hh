//----------------------------------*-C++-*----------------------------------//
// Release.hh
// Randy M. Roberts
// Sun Aug 22 08:26:12 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function for POOMA_MT library
//---------------------------------------------------------------------------//

#ifndef __POOMA_MT_Release_hh__
#define __POOMA_MT_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of POOMA_MT; 
// this can be used to get exact version information in codes that 
// use POOMA_MT
// 
//===========================================================================//

#include <string>
// #include <POOMA_MT/config.h>

namespace rtt_POOMA_MT 
{
    const std::string release();
}

#endif                          // __POOMA_MT_Release_hh__

//---------------------------------------------------------------------------//
//                              end of POOMA_MT/Release.hh
//---------------------------------------------------------------------------//
