//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file  ds++/SystemCall.cc
 * \brief Implementation for the Draco wrapper for system calls. This routine
 *        attempts to hide differences between Unix/Windows system calls.
 * \note  Copyright (C) 2016-2020 Triad National Security, LLC.
 *        All rights reserved. */
//----------------------------------------------------------------------------//

#include "SystemCall.hh"
#include "Assert.hh"

#include <cerrno>  // errno
#include <cstdio>  // remove()
#include <cstdlib> // _fullpath
#include <cstring> // strncpy()
#if defined UNIX || defined MINGW
#include <sys/param.h> // MAXPATHLEN
#include <unistd.h>    // gethostname
#endif
#ifdef WIN32
#include <direct.h>  // _getcwd
#include <process.h> // _getpid
#include <sstream>
#endif

#include <iostream>

namespace rtt_dsxx {

//============================================================================//
// FREE FUNCTIONS
//============================================================================//

//----------------------------------------------------------------------------//
/*! \brief Wrapper for system dependent hostname call.
 *
 * Windows:
 *     \c HOST_NAME_MAX set to \c MAX_COMPUTERNAME_LENGTH in config.h
 *
 * Catamount systems:
 *     \c HOST_NAME_MAX hard coded by CMake in config.h
 *
 * Unix/Linux:
 *     \c HOST_NAME_MAX loaded from \<climit\>
 *
 * Mac OSX:
 *     \c HOST_NAME_MAX set to \c _POSIX_HOST_NAME_MAX in config.h
 */
std::string draco_gethostname() {
// Windows: gethostname from <winsock2.h>
#ifdef WIN32
  char hostname[HOST_NAME_MAX];
  int err = gethostname(hostname, sizeof(hostname));
  if (err)
    strncpy(hostname, "gethostname() failed", HOST_NAME_MAX);
  return std::string(hostname);

#else

// Linux: gethostname from <unistd.h>
#ifdef HAVE_GETHOSTNAME
  char hostname[HOST_NAME_MAX];
  int err = gethostname(hostname, HOST_NAME_MAX);
  if (err)
    strncpy(hostname, "gethostname() failed", HOST_NAME_MAX);
  return std::string(hostname);

// Catamount systems do not have gethostname().
#else
  return std::string("Host (unknown)");
#endif
#endif
} // draco_gethostname

//----------------------------------------------------------------------------//
/*! \brief Wrapper for system dependent pid call..
 *
 * Catamount systems do not have getpid().  This function will return -1.
 */
int draco_getpid() {
#ifdef WIN32
  int i = _getpid();
  return i;
#else

#ifdef HAVE_GETHOSTNAME
  return getpid();
#else
  // Catamount systems do not have getpid().  This function will return -1.
  return -1;
#endif
#endif
} // draco_pid

//----------------------------------------------------------------------------//
/*! \brief Wrapper for system dependent getcwd call.
 *
 *  This should always return a trailing directory separator.
 */
std::string draco_getcwd() {
// Identify the current working directory.
#ifdef WIN32
  char *buffer;
  Insist((buffer = _getcwd(nullptr, 0)) != nullptr,
         std::string("getcwd failed: " + std::string(strerror(errno))));
  std::string cwd(buffer, buffer + strnlen(buffer, MAXPATHLEN));
  free(buffer);
#else
  char curr_path[MAXPATHLEN];
  curr_path[0] = '\0';
  Insist(getcwd(curr_path, MAXPATHLEN) != nullptr,
         std::string("getcwd failed: " + std::string(strerror(errno))));
  std::string cwd(curr_path);
#endif

  // Ensure that the last character is a directory separator.
  if (cwd.at(cwd.length() - 1) != rtt_dsxx::dirSep)
    cwd = cwd + rtt_dsxx::dirSep;

  return cwd;
}

//----------------------------------------------------------------------------//
/*! \brief Wrapper for system dependent stat call.
 *
 * Windows
 *    http://msdn.microsoft.com/en-us/library/14h5k7ff%28v=vs.80%29.aspx
 * Linux
 *    http://en.wikipedia.org/wiki/Stat_%28system_call%29
 */
draco_getstat::draco_getstat(std::string const &fqName)
    : stat_return_code(0), buf() {
#ifdef WIN32
  filefound = true;
  /*! \note If path contains the location of a directory, it cannot
     * contain a trailing backslash. If it does, -1 will be returned and
     * errno will be set to ENOENT. */
  std::string clean_fqName;
  if (fqName[fqName.size() - 1] == rtt_dsxx::WinDirSep ||
      fqName[fqName.size() - 1] == rtt_dsxx::UnixDirSep)
    clean_fqName = fqName.substr(0, fqName.size() - 1);
  else
    clean_fqName = fqName;

  stat_return_code = _stat(clean_fqName.c_str(), &buf);
  if (stat_return_code != 0) {
    switch (errno) {
    case ENOENT: // file not found
      filefound = false;
      break;
    case EINVAL:
      Insist(stat_return_code == 0, "invalid parameter given to _stat.");
      break;
    default:
      /* should never get here */
      Insist(stat_return_code == 0, "_stat returned an error.");
    }
  }

  if (filefound) {
    // Handle to directory
    HANDLE hFile = ::FindFirstFile(clean_fqName.c_str(), &FileInformation);
    // sanity check
    Insist(hFile != INVALID_HANDLE_VALUE, "Invalid file handle.");
    // Close handle
    ::FindClose(hFile);
  }

#else
  stat_return_code = stat(fqName.c_str(), &buf);
#endif
}

//----------------------------------------------------------------------------//
//! Is this a regular file?
bool draco_getstat::isreg() {
#ifdef WIN32
  bool b = FileInformation.dwFileAttributes & FILE_ATTRIBUTE_NORMAL;
  return filefound && b;
#else
  bool b = S_ISREG(buf.st_mode);
  return b;
#endif
}

//----------------------------------------------------------------------------//
//! Is this a directory?
bool draco_getstat::isdir() {
#ifdef WIN32
  bool b = FileInformation.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY;
  return filefound && b;
#else
  bool b = S_ISDIR(buf.st_mode);
  return b;
#endif
}

//----------------------------------------------------------------------------//
//! Is a Unix permission bit set?
#ifdef WIN32
bool draco_getstat::has_permission_bit(int /*mask*/) {
  Insist(isreg(), "Can only check permission bit for regular files.");
  Insist(false,
         "draco_getstat::hsa_permission_bit() not implemented for WIN32");
  return false;
}
#else
bool draco_getstat::has_permission_bit(int mask) {
  Insist(isreg(), "Can only check permission bit for regular files.");
  // check execute bit (buf.st_mode & 0111)
  return (buf.st_mode & mask);
}
#endif

//----------------------------------------------------------------------------//
/*!
 * \brief Wrapper for system dependent realpath call.
 */
std::string draco_getrealpath(std::string const &path) {
  char buffer[MAXPATHLEN]; // _MAX_PATH
#ifdef WIN32
  // http://msdn.microsoft.com/en-us/library/506720ff%28v=vs.100%29.aspx
  Insist(_fullpath(buffer, path.c_str(), MAXPATHLEN) != nullptr,
         "Invalid path.");
  std::string retVal(buffer);
#else
  Insist((realpath(path.c_str(), buffer)) != nullptr, "Invalid path.");
  // realpath trims the trailing slash, append now.
  std::string retVal(buffer);
  retVal += std::string(&dirSep, 1);
#endif
  return retVal;
}

//----------------------------------------------------------------------------//
/*!
 * \brief Make a directory
 * boost::filesystem::create_directories("/tmp/a/b/c");
 * \todo Do we need a permissions argument?
 */
void draco_mkdir(std::string const &path) {
#ifdef WIN32

  draco_getstat dirInfo(path);
  if (!dirInfo.isdir()) {
    /*! \note If path contains the location of a directory, it cannot contain a
     * trailing backslash. If it does, -1 will be returned and errno will be set
     * to ENOENT. */
    std::string clean_fqName;
    if (path[path.size() - 1] == rtt_dsxx::WinDirSep ||
        path[path.size() - 1] == rtt_dsxx::UnixDirSep)
      clean_fqName = path.substr(0, path.size() - 1);
    else
      clean_fqName = path;

    // make the directory
    int errorCode = _mkdir(clean_fqName.c_str());
    if (errorCode == -1) {
      switch (errno) {
      case EEXIST:
        Insist(errno != EEXIST, "ERROR: Unable to create directory," +
                                    clean_fqName +
                                    ", because it already exists.");
        break;
      case ENOENT:
        Insist(errno != ENOENT, "ERROR: Unable to create directory," +
                                    clean_fqName +
                                    ", because the path is not found.");
        break;
      default:
        /* should never get here */
        Insist(errno == 0,
               "ERROR: Unable to create directory, " + clean_fqName);
      }
    }
  }
#else
  mkdir(path.c_str(), 0700);
// S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
#endif
}

//----------------------------------------------------------------------------//
/*!
 * \brief Delete a single file or directory (not recursive)
 *
 * \sa wdtOpRemove
 *
 * Sample implementation for Win32 (uses Win API which I don't want to do)
 * http://forums.codeguru.com/showthread.php?239271-Windows-SDK-File-System-How-to-delete-a-directory-and-subdirectories
 * http://stackoverflow.com/questions/1468774/recursive-directory-deletion-with-win32
 * http://msdn.microsoft.com/en-us/library/aa365488%28VS.85%29.aspx
 *
 * Sample implementation for Unix
 * http://www.linuxquestions.org/questions/programming-9/deleting-a-directory-using-c-in-linux-248696/
 *
 * Consider using Boost.FileSystem
 */
void draco_remove(std::string const &dirpath) {
  // remove() works for all unix items but only for files
  // (not directories) for windows.
  if (draco_getstat(dirpath).isdir()) {
#ifdef WIN32
    // Clear any special directory attributes.
    bool ok = ::SetFileAttributes(dirpath.c_str(), FILE_ATTRIBUTE_NORMAL);
    if (!ok) {
      int myerr = ::GetLastError();
      std::ostringstream msg;
      msg << "ERROR: File attribute not normal. myerr = " << myerr
          << ", file = " << dirpath;
      Insist(ok, msg.str());
    }

    // Delete directory
    ok = ::RemoveDirectory(dirpath.c_str());
    if (!ok) {
      int myerr = ::GetLastError();
      std::ostringstream msg;
      msg << "ERROR: Error deteting file, myerr = " << myerr
          << ", file = " << dirpath;
      Insist(ok, msg.str());
    }
#else
    remove(dirpath.c_str());
#endif
  } else {
    remove(dirpath.c_str());
  }
  // If the file still exists this check will fail.
  Ensure(!draco_getstat(dirpath).valid());
  return;
}

} // end of namespace rtt_dsxx

//----------------------------------------------------------------------------//
// end of SystemCall.cc
//----------------------------------------------------------------------------//
