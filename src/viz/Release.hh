//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/Release.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 21 16:29:46 2000
 * \brief  Release function for the viz library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __viz_Release_hh__
#define __viz_Release_hh__

//===========================================================================//
/*!
 * \page viz_overview Overview of the viz package
 * \version 1_0_0
 * \author T.M. Evans and J. McGhee
 * 
 * The viz package contains translators to visualization software.  Currently
 * supported are output translators for EnSight.
 *
 * More translators will come on line as need arises.
 */
//===========================================================================//
/*!
 * \namespace rtt_viz
 *
 * \brief Namespace that contains the viz package classes and variables.
 *
 */
//===========================================================================//

#include <string>

namespace rtt_viz 
{
    const std::string release();
}

#endif                          // __viz_Release_hh__

//---------------------------------------------------------------------------//
//                           end of viz/Release.hh
//---------------------------------------------------------------------------//
