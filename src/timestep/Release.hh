//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   timestep/Release.hh
 * \author Thomas M. Evans
 * \date   Mon Apr 19 21:35:59 2004
 * \brief  Release function for the timestep library
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef timestep_Release_hh
#define timestep_Release_hh

#include <string>

//===========================================================================//
/*!
 * \namespace rtt_timestep
 *
 * \brief Namespace that contains the timestep package classes and variables.
 *
 * Provides namespace protection for the Draco (RTT) time step control
 * utilities.
 *
 * \sa The ts_manager and ts_advisor classes provide most of the
 * functionality of the namespace. The <a href="./index.html">Main Page</a>
 * presents a summary of the capabilities provided within the namespace.
 */
/*!
 * \example timestep/test/main.cc
 *
 * The following code provides an example of how to use the timestep manager
 * utility. It includes a dummy_package for use with the manager. Also
 * included isis a sample output from the test_timestep example.  It contains
 * representative output from most of the printing and summary I/O utilities.
 *
 * \include timestep/test/test_timestep.hh
 * \include timestep/test/test_timestep.cc
 * \include timestep/test/dummy_package.hh
 * \include timestep/test/test_timestep_pt.cc
 * \include timestep/test/dummy_package.cc
 * \include timestep/test/test_timestep.out
 */
//===========================================================================//

namespace rtt_timestep
{
    const std::string release();
}

#endif // timestep_Release_hh

//---------------------------------------------------------------------------//
//                        end of timestep/Release.hh
//---------------------------------------------------------------------------//
