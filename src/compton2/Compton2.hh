//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   compton2/Compton2.hh
 * \author Andrew Till
 * \date   11 May 2020
 * \brief  Header file for compton interface
 * \note   Copyright (C) 2017-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef rtt_compton2_Compton2_hh
#define rtt_compton2_Compton2_hh

// C++ standard library dependencies
#include <iostream>
#include <memory>
#include <vector>

namespace rtt_compton2 {
//============================================================================//
/*!
 * \class Compton2
 *
 * \brief Provides access to relativistic Compton scattering angle and
 *        multigroup frequency distributions from csk data files
 *
 * This interface class allows the client to:
 * 1) access (interpolate) data from existing multigroup CSK_generator libraries
 * 2) build new multigroup libraries from existing CSK_generator pointwise
      libraries
 * 3) obtain auxiliary information for existing multigroup libraries
 *    (electron temperature bounds, frequency group structures, etc)
 *
 */

/*!
 * \example compton2/test/tCompton2.cc
 *
 * This unit test demonstrates the two methods for constructing a Compton
 * object, and exercises all routines for interpolation and data access.
*/

class Compton2 {};

} // namespace rtt_compton2

#endif // rtt_compton2_Compton2_hh

//----------------------------------------------------------------------------//
// End compton2/Compton2.hh
//----------------------------------------------------------------------------//
