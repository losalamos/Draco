//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/DACS_Process.cc
 * \author Gabriel M. Rockefeller
 * \date   Thu Jul 21 13:09:19 2011
 * \brief  DACS_Process implementation.
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DACS_Process.hh"

#include "device/config.h"
#include "ds++/Assert.hh"

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <sys/param.h>
#include <sys/stat.h>

namespace rtt_device
{

static std::string canonical_fname(std::string const & path);
static std::string verbose_error(std::string const & message);

//---------------------------------------------------------------------------//
/*! \brief DACS_Process constructor.
 *
 * DACS_Process saves the full filename of the accel-side binary.  The
 * constructor tries to interpret the specified filename in three ways: as an
 * absolute path, as a path relative to the current working directory, and as
 * a path relative to a default location (default_ppe_bindir).  If the binary
 * can't be found in any of those locations, the constructor throws an
 * exception.
 */
DACS_Process::DACS_Process(const std::string & filename)
{
    // Identify the current working directory.
    char curr_path[MAXPATHLEN];
    Insist(getcwd(curr_path, MAXPATHLEN) != NULL,
           verbose_error("getcwd failed: " + std::string(strerror(errno))));

    const std::string cwd(curr_path);

    const std::string fullname1 = cwd + "/" + filename;
    const std::string fullname2 = default_ppe_bindir + "/" + filename;

    // Interpret filename as an absolute pathname...
    std::string fullname = canonical_fname(filename);
    if (fullname.empty())
    {
        // ... as a pathname relative to the current directory...
        fullname = canonical_fname(fullname1);
        if (fullname.empty())
        {
            // ... and as a pathname relative to default_ppe_bindir.
            fullname = canonical_fname(fullname2);
            Insist(!fullname.empty(),
                   verbose_error("Couldn't stat " + filename +
                                 " or " + fullname1 +
                                 " or " + fullname2 + ": " +
                                 std::string(strerror(errno))));
        }
    }

    accel_filename = fullname;
} // DACS_Process::DACS_Process

//---------------------------------------------------------------------------//
/*! \brief Launch the accel-side process.
 *
 * start invokes dacs_de_start to launch the accel-side binary.
 */
void DACS_Process::start(const de_id_t & de_id)
{
    void * vf =
        reinterpret_cast<void*>(const_cast<char *>(accel_filename.c_str()));

    DACS_ERR_T err = dacs_de_start(de_id,
                                   vf,
                                   NULL,
                                   NULL,
                                   DACS_PROC_LOCAL_FILE,
                                   &pid);

    Insist(err == DACS_SUCCESS,
           verbose_error(accel_filename + ": " +
                         std::string(dacs_strerror(err))));
} // DACS_Process::start

//---------------------------------------------------------------------------//
/*! \brief Return the canonical form of a pathname.
 *
 * DACS has trouble with pathnames containing symlinks.  canonical_fname
 * checks to make sure that path exists (via stat) and then invokes realpath
 * to expand symbolic links.  If stat returns an error, canonical_fname
 * returns an empty string (and the DACS_Device constructor will try another
 * location).
 */
std::string canonical_fname(std::string const &path)
{
    struct stat buf;
    int err = stat(path.c_str(), &buf);
    if (err) return std::string();
    
    char npath[PATH_MAX]; npath[0] = '\0';
    Insist(realpath(path.c_str(), npath) != NULL,
           verbose_error("realpath failed on " + path));

    return std::string(npath);
} // canonical_fname

//---------------------------------------------------------------------------//
/*! \brief Add hostname and pid to error messages.
 *
 * Several of the errors that might be reported by DACS_Device could be
 * specific to one or a few nodes (filesystems not mounted, etc.).
 * verbose_error adds the hostname and pid to error messages.
 */
std::string verbose_error(std::string const & message)
{
    char hostname[HOST_NAME_MAX];
    int err = gethostname(hostname, HOST_NAME_MAX);
    if (err) strncpy(hostname, "gethostname() failed", HOST_NAME_MAX);

    std::ostringstream errstr;
    errstr << "Host " << hostname << ", PID " << getpid() << ": "
           << message;

    return errstr.str();
} // verbose_error

} // end namespace rtt_device

//---------------------------------------------------------------------------//
//                 end of DACS_Process.cc
//---------------------------------------------------------------------------//
