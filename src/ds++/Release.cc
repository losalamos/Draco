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

#include "Release.hh"
#include "ds++/config.h"
#include <sstream>
#include <cstring> // memcpy

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

//! This version can be called by Fortran and wraps the C++ version.
extern "C" void ec_release( char * release_string, size_t maxlen )
{
    std::string tmp_rel = rtt_dsxx::release();
    if( tmp_rel.size() >= static_cast<size_t>(maxlen) )
    {
        tmp_rel = tmp_rel.substr(0,maxlen-1);
    }
    std::memcpy(release_string,tmp_rel.c_str(),tmp_rel.size()+1);
    return;    
}

//---------------------------------------------------------------------------//
//                              end of Release.cc
//---------------------------------------------------------------------------//
