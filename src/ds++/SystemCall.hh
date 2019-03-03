//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/SystemCall.hh
 * \brief  Wrapper for system calls. Hide differences between Unix/Windows
 *         system calls.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef rtt_dsxx_SystemCall_hh
#define rtt_dsxx_SystemCall_hh

#include "ds++/config.h"
#include <string>
#ifdef WIN32
#include <WinSock2.h>
#include <Windows.h>
#include <sys/types.h>
#endif
#include <sys/stat.h> // stat (UNIX) or _stat (WIN32)

namespace rtt_dsxx {

//! Character used as path separator.
char const WinDirSep = '\\';
char const UnixDirSep = '/';
#ifdef _MSC_VER
char const dirSep = WinDirSep;
std::string const exeExtension(".exe");
#else
char const dirSep = UnixDirSep;
std::string const exeExtension("");
#endif

//===========================================================================//
// General discussion.  See .cc file for detailed implementation discussion
// (mostly Linux vs. Windows issues).
//===========================================================================//

/*! \section HOST_NAME_MAX HOST_NAME_MAX
 *
 * The selection of a value for HOST_NAME_MAX is completed by
 * ds++/CMakeLists.txt and ds++/config.h.in.
 *
 * - For most Linux platforms, \c HOST_NAME_MAX is defined in \c \<limits.h\>.
 *   However, according to the POSIX standard, \c HOST_NAME_MAX is a
 *   \em possibly \em indeterminate definition meaning that it
 *
 * \note ...shall be omitted from \c \<limits.h\> on specific implementations
 *       where the corresponding value is equal to or greater than the stated
 *       minimum, but is unspecified.
 *
 * - The minumum POSIX guarantee is \c HOST_NAME_MAX = \c 256.
 * - An alternate value used by some Unix systems is \c MAXHOSTNAMELEN as
 *   defined in \c \<sys/param.h\>
 * - On Windows, the variable \c MAX_COMPUTERNAME_LENGTH from \c \<windows.h\>
 *   can be used. See http://msdn.microsoft.com/en-us/library/windows/desktop/ms738527%28v=vs.85%29.aspx
 *  - On Mac OSX, we use \c _POSIX_HOST_NAME_MAX.
 */

//===========================================================================//
// FREE FUNCTIONS
//===========================================================================//

//! Return the local hostname
DLL_PUBLIC_dsxx std::string draco_gethostname(void);

//! Return the local process id
DLL_PUBLIC_dsxx int draco_getpid(void);

//! Return the current working directory
DLL_PUBLIC_dsxx std::string draco_getcwd(void);

//! Return the stat value for a file
class DLL_PUBLIC_dsxx draco_getstat {
private:
  int stat_return_code;
#ifdef WIN32
  struct _stat buf;
  bool filefound;
  WIN32_FIND_DATA FileInformation; // Additional file information
#else
  struct stat buf;
#endif

public:
  //! constructor
  explicit draco_getstat(std::string const &fqName);
  //! If the call to stat failed, this function will return false.
  bool valid(void) { return stat_return_code == 0; };
  bool isreg(void);
  bool isdir(void);
  int errorCode(void) { return stat_return_code; }
  /*!
   * \brief Determine if the file has the requested permission bits set.
   * \note The leading zero for the mask is important.
   */
  bool has_permission_bit(int mask = 0777);
};

//! Use Linux realpath to resolve symlinks
DLL_PUBLIC_dsxx std::string draco_getrealpath(std::string const &path);

//! Create a directory
DLL_PUBLIC_dsxx void draco_mkdir(std::string const &path);

/*!
 * \brief Remove file or directory (not recursive)
 *
 * For recursive directory delete, see path.hh's walk_directory_tree and
 * the functor wdtOpRemove.
 */
DLL_PUBLIC_dsxx void draco_remove(std::string const &path);

} // namespace rtt_dsxx

#endif // rtt_dsxx_SystemCall_hh

//---------------------------------------------------------------------------//
// end of SystemCall.hh
//---------------------------------------------------------------------------//
