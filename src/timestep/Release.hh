//----------------------------------*-C++-*----------------------------------//
// Release.hh
// Randy M. Roberts
// Fri Aug 27 10:33:26 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function for timestep library
//---------------------------------------------------------------------------//

#ifndef __timestep_Release_hh__
#define __timestep_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of timestep; 
// this can be used to get exact version information in codes that 
// use timestep
// 
//===========================================================================//

#include <string>

namespace rtt_timestep 
{
    const std::string release();
}

#endif                          // __timestep_Release_hh__

//---------------------------------------------------------------------------//
//                              end of timestep/Release.hh
//---------------------------------------------------------------------------//
