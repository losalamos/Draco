//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   stdheaders/Release.hh
 * \author Michelle L. Murillo
 * \date   Wed Jun 21 17:46:55 2000
 * \brief  Release function for the stdheaders library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __stdheaders_Release_hh__
#define __stdheaders_Release_hh__

//===========================================================================//
/*!
 * \mainpage Overview of the stdheaders package
 * \version 1_0_0
 * \author Michelle Murillo and Tom Evans
 * 
 * stdheaders provides C++ Standard defined headers for the ANSI C-STD
 * library.  The C++ standard stipulates that c-std libraries, ie. math.h,
 * should sit in the std namespace and be called using c(library name).
 * Thus, cos from math.h is called using std::cos from cmath.
 *
 * stdheaders is turned on by the --enable-draco-stdhdrs option.
 * 
 */
//===========================================================================//
/*!
 * \namespace rtt_stdheaders
 *
 * \brief Namespace that contains the stdheaders package classes and variables.
 *
 */
//===========================================================================//

#include <string>

namespace rtt_stdheaders 
{
    const std::string release();
}

#endif                          // __stdheaders_Release_hh__

//---------------------------------------------------------------------------//
//                           end of stdheaders/Release.hh
//---------------------------------------------------------------------------//
