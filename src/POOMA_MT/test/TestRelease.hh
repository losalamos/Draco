//----------------------------------*-C++-*----------------------------------//
// Release.hh
// Randy M. Roberts
// Wed Aug 25 16:27:15 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function for POOMA_MT/test programs.
//---------------------------------------------------------------------------//

#ifndef __POOMA_MT_test_Release_hh__
#define __POOMA_MT_test_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of POOMA_MT/test; 
// this can be used to get exact version information in codes that 
// use POOMA_MT
// 
//===========================================================================//

#include <string>

namespace rtt_POOMA_MT_test 
{
    const std::string release();
}

#endif                          // __POOMA_MT_test_Release_hh__

//---------------------------------------------------------------------------//
//                              end of POOMA_MT/test/Release.hh
//---------------------------------------------------------------------------//
