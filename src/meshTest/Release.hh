//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/Release.hh
 * \author Randy M. Roberts
 * \date   Wed Sep 22 15:24:01 1999
 * \brief  Release function for meshTest system
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_Release_hh__
#define __meshTest_Release_hh__

//===========================================================================//
/*!
 * \page meshTest_overview Overview of the meshTest system
 * \version 1_0_1
 * \author Shawn Pautz, Randy M. Roberts
 *
 * The meshTest system is used to test a Mesh MT's services.
 */
//===========================================================================///
/*!
 * \namespace rtt_meshTest
 *
 * \brief Namespace that contains the meshTest package classes and variables.
 *
 * The rtt_meshTest namespace contains all classes and variables necessary to
 * use the meshTest system.  This system is used to test Mesh MT's adherence
 * to the MT service requirements.
 */
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
