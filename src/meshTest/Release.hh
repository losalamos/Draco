//----------------------------------*-C++-*----------------------------------//
// Release.hh
// Randy M. Roberts
// Fri Aug 20 13:04:58 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function for meshTest library
//---------------------------------------------------------------------------//

#ifndef __meshTest_Release_hh__
#define __meshTest_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of meshTest; 
// this can be used to get exact version information in codes that 
// use meshTest
// 
//===========================================================================//

#include <string>

namespace rtt_meshTest 
{
    const std::string release();
}

#endif                          // __meshTest_Release_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/Release.hh
//---------------------------------------------------------------------------//
