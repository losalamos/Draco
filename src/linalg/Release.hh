//----------------------------------*-C++-*----------------------------------//
// Release.hh
// Randy M. Roberts
// Thu Jul 15 09:44:12 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function for linalg library
//---------------------------------------------------------------------------//

#ifndef __linalg_Release_hh__
#define __linalg_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of linalg; 
// this can be used to get exact version information in codes that 
// use linalg
// 
//===========================================================================//

#include <string>

namespace rtt_linalg 
{
    const std::string release();
}

#endif                          // __linalg_Release_hh__

//---------------------------------------------------------------------------//
//                              end of linalg/Release.hh
//---------------------------------------------------------------------------//
