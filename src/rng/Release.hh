//----------------------------------*-C++-*----------------------------------//
// Release.hh
// Thomas M. Evans
// Thu May 27 15:24:01 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function for rng library
//---------------------------------------------------------------------------//

#ifndef __rng_Release_hh__
#define __rng_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of rng; 
// this can be used to get exact version information in codes that 
// use rng
// 
//===========================================================================//

#include <string>

namespace rtt_rng 
{
    const std::string release();
}

#endif                          // __rng_Release_hh__

//---------------------------------------------------------------------------//
//                              end of rng/Release.hh
//---------------------------------------------------------------------------//
