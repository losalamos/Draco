//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   units/Release.hh
 * \author Randy M. Roberts
 * \date   Wed Feb  9 10:57:46 2000
 * \brief  Release function for the units library
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __units_Release_hh__
#define __units_Release_hh__

//===========================================================================//
/*!
 * \page units_overview Overview of the units package
 * \version 2.0
 * \author Kelly Thompson
 *
 * This package provides a mechanism to keep track of, compare and convert
 * between unit systems.  It also provides access to physical constants
 * (speed of light, pi, avogadro's number, etc.)  The general usage is:
 *
 * UnitSystem SI_Units( UnitSystemType:SI() );
 * UnitSystem X4_Units( UnitSystemType:X4() );
 * UnitSystem my_Units( UnitSystemType().L( rtt_units::L_m )
 *                                      .M( rtt_units::M_g )
 *                                      .t( rtt_units::t_us ) );
 * PhysicalConstants pc( X4_Units );
 *
 * x4_speed_of_light = pc.c(); // cm/sh
 *
 * my_energy = 5.0; // g * m^2 / us^2
 *
 * my_energy_si = my_energy * SI_Units.e() / my_Units.e(); // J = kg * m^2 / s^2
 * 
 */
//===========================================================================//
/*!
 * \namespace rtt_units
 *
 * \brief Namespace that contains the units package classes and variables.
 *
 */
//===========================================================================//

#include <string>

namespace rtt_units 
{
    const std::string release();
}

#endif                          // __units_Release_hh__

//---------------------------------------------------------------------------//
//                           end of units/Release.hh
//---------------------------------------------------------------------------//
