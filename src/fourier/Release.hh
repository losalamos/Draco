//----------------------------------*-C++-*----------------------------------//
// Release.hh
// John Gulick
// Tue Aug 24 13:08:51 1999
// $Id$
//---------------------------------------------------------------------------//
// @> Release function for fourier library
//---------------------------------------------------------------------------//

#ifndef __fourier_Release_hh__
#define __fourier_Release_hh__

//===========================================================================//
// namespace version - 
//
// Purpose : Return the version of fourier; 
// this can be used to get exact version information in codes that 
// use fourier
// 
//===========================================================================//

#include <string>

namespace rtt_fourier 
{
    const std::string release();
}

#endif                          // __fourier_Release_hh__

//---------------------------------------------------------------------------//
//                              end of fourier/Release.hh
//---------------------------------------------------------------------------//
