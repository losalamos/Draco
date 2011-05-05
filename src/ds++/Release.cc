//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ../Release.cc
 * \author Thomas Evans
 * \date   Thu Jul 15 09:31:44 1999
 * \brief  Provides the function definition for Release.
 * \note   Copyright (C) 1999-2010 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <sstream>
#include "Release.hh"
#include "ds++/config.h"

namespace rtt_dsxx
{

// function definition for Release, define the local version number for
// this library in the form ds_#.#.# in pkg_version variable
const std::string release()
{
    std::ostringstream pkg_release;
    pkg_release << "Draco-"
                << DRACO_VERSION_MAJOR << "_"
                << DRACO_VERSION_MINOR << "_"
                << DRACO_VERSION_PATCH ;
    return pkg_release.str();
}

}  // end of rtt_dsxx

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
