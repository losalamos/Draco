//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   ds++/SystemCall.hh
 * \brief Wrapper for system calls. Hide differences between Unix/Windows system calls.
 * \note Copyright (C) 2016-2020 Triad National Security, LLC.  All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef rtt_dsxx_SystemCall_hh
#define rtt_dsxx_SystemCall_hh

#include "ds++/config.h"
#include <string>
#ifdef WIN32
#define _WINSOCKAPI_
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

//================================================================================================//
// FREE FUNCTIONS
//================================================================================================//

//! Return the local hostname
std::string draco_gethostname();

//! Return the local process id
int draco_getpid();

//! Return the current working directory
std::string draco_getcwd();

//! Return the stat value for a file
class draco_getstat {
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
  bool valid() { return stat_return_code == 0; };
  bool isreg();
  bool isdir();
  int errorCode() { return stat_return_code; }
  /*!
   * \brief Determine if the file has the requested permission bits set.
   * \note The leading zero for the mask is important.
   */
  bool has_permission_bit(int mask = 0777);
};

//! Use Linux realpath to resolve symlinks
std::string draco_getrealpath(std::string const &path);

//! Create a directory
void draco_mkdir(std::string const &path);

/*!
 * \brief Remove file or directory (not recursive)
 *
 * For recursive directory delete, see path.hh's walk_directory_tree and the functor wdtOpRemove.
 */
void draco_remove(std::string const &path);

} // namespace rtt_dsxx

/*!
 * \brief set and unset environment variables.
 */
#ifdef MSVC
#define draco_unsetenv(k) _putenv_s(k, "")
#define draco_setenv(k, v) _putenv_s(k, v)
#else
#define draco_unsetenv(k) unsetenv(k)
#define draco_setenv(k, v) setenv(k, v, 1)
#endif

#endif // rtt_dsxx_SystemCall_hh

//------------------------------------------------------------------------------------------------//
// end of SystemCall.hh
//------------------------------------------------------------------------------------------------//
