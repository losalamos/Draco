//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/path.hh
 * \author Kelly Thompson
 * \date   Thu Sep 27 8:49:19 2010
 * \brief  Encapsulate path information (path separator, etc.)
 * \note   Copyright © 2010 Los Alamos National Security, LLC
 *         All rights reserved.
 * \version $Id$
 */
//---------------------------------------------------------------------------//

/*!
 * \bug Consider replacing path.cc and path.hh with Boost FileSystem.
 */

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

//---------------------------------------------------------------------------//
/*! \brief Report the current working directory
 */
std::string currentPath(void);

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
std::string getFilenameComponent( std::string const & fqName,
                                  FilenameComponent fc );
    
} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of path.hh
//---------------------------------------------------------------------------//
