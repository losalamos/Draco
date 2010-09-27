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

namespace rtt_dsxx
{

//! Character used as path separator.
char const WinDirSep  = '\\';
char const UnixDirSep = '/';
#ifdef _MSC_VER
char const dirSep = WinDirSep;
#else
char const dirSep = UnixDirSep;
#endif

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of path.hh
//---------------------------------------------------------------------------//
