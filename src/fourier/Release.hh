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
/*!
 * \page overview Overview of the Fourier Analysis Package
 *
 * This package is used to do Fourier analysis.  The primary classes are
 * rtt_fourier::fourier, rtt_fourier::matrix, rtt_fourier::gauss.  It uses
 * the mtl library to do eigensolves.
 *
 * This is basically how its used.  Register a class that has the creates a
 * matrix of type X.  This ....  
*/
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
