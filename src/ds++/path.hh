//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/path.hh
 * \brief  Encapsulate path information (path separator, etc.)
 * \note   Copyright © 2010 Los Alamos National Security, LLC
 *         All rights reserved.
 * \version $Id$
 */
//---------------------------------------------------------------------------//

#ifndef __dsxx_path_hh__
#define __dsxx_path_hh__

/*!
 * \bug Consider replacing path.cc and path.hh with Boost FileSystem.
 */

#include <ds++/config.h>
#include <string>

namespace rtt_dsxx
{

//! Character used as path separator.
char const WinDirSep  = '\\';
char const UnixDirSep = '/';
#ifdef _MSC_VER
char const dirSep = WinDirSep;
std::string const exeExtension( ".exe" );
#else
char const dirSep = UnixDirSep;
std::string const exeExtension( "" );
#endif

enum FilenameComponent
{
    FC_PATH,       //!< Extract path portion of fully qualified filename
    FC_ABSOLUTE,
    FC_NAME,       //!< Extract filename portion (w/o path).
    FC_EXT,
    FC_NAME_WE,
    FC_REALPATH,
    FC_LASTVALUE
};

//---------------------------------------------------------------------------//
/*! \brief Get a specific component of a full filename.
 */
DLL_PUBLIC std::string getFilenameComponent( std::string const & fqName,
                                             FilenameComponent   fc );

//---------------------------------------------------------------------------//
/*! \brief Does the file exist?
 */
DLL_PUBLIC bool fileExists( std::string const & filename );

} // end namespace rtt_dsxx

#endif // __dsxx_path_hh__

//---------------------------------------------------------------------------//
// end of path.hh
//---------------------------------------------------------------------------//
