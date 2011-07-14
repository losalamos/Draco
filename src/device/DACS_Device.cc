//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/DACS_Device.cc
 * \author Gabriel M. Rockefeller
 * \date   Thu Jun 16 14:55:48 2011
 * \brief  DACS_Device implementation.
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DACS_Device.hh"

#include "device/config.h"
#include "ds++/Assert.hh"

#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include <string>
#include <sys/param.h>
#include <sys/stat.h>
#include <unistd.h>

namespace rtt_device
{

std::string DACS_Device::accel_filename;

std::string canonical_fname(std::string const & path);
std::string verbose_error(std::string const & message);

//---------------------------------------------------------------------------//
/*! \brief DACS_Device default constructor.
 *
 * DACS_Device is a singleton, so the constructor is called at most once, the
 * first time any other code invokes DACS_Device::instance.  The constructor
 * initializes DACS (via dacs_init), reserves a Cell (via
 * dacs_reserve_children) and launches the accel-side binary (via
 * dacs_de_start).  If no accel-side binary has been specified (via
 * DACS_Device::init), the constructor throws an exception.
 */
DACS_Device::DACS_Device() : de_id(0), pid(0)
{
    Insist(!accel_filename.empty(), "No accel filename specified");

    // Initialize DACS.
    DACS_ERR_T err = dacs_init(DACS_INIT_FLAGS_NONE);
    Insist(err == DACS_SUCCESS, verbose_error(dacs_strerror(err)));

    // Get the number of available Cell BEs.
    uint32_t num_children(0);
    err = dacs_get_num_avail_children(DACS_DE_CBE, &num_children);
    Insist(err == DACS_SUCCESS, verbose_error(dacs_strerror(err)));

    Insist(num_children > 0, "No child CBE available");

    // Reserve one Cell BE.
    uint32_t num_reserve(1);
    err = dacs_reserve_children(DACS_DE_CBE, &num_reserve, &de_id);
    Insist(err == DACS_SUCCESS, verbose_error(dacs_strerror(err)));

    Insist(num_reserve == 1, "Failed to reserve children");

    // Start the accel-side process.
    void * vf =
        reinterpret_cast<void*>(const_cast<char *>(accel_filename.c_str()));
    err = dacs_de_start(de_id, vf, NULL, NULL, DACS_PROC_LOCAL_FILE, &pid);
    Insist(err == DACS_SUCCESS,
           verbose_error(accel_filename + ": " +
                         std::string(dacs_strerror(err))));

} // DACS_Device::DACS_Device

//---------------------------------------------------------------------------//
/*! \brief DACS_Device destructor.
 *
 * DACS_Device is a singleton; the destructor is called at program
 * termination.  The destructor waits for the accel-side process to complete
 * and then calls dacs_exit.  Throwing an exception inside the destructor
 * would interfere with stack unwinding already in progress (from, for
 * example, an earlier exception), so errors are reported to stderr.
 */
DACS_Device::~DACS_Device()
{
    // Wait for the accel-side process to exit.
    int32_t exit_status;
    DACS_ERR_T err = dacs_de_wait(de_id,
                                  pid,
                                  &exit_status);
    if (err != DACS_STS_PROC_FINISHED)
        std::cerr << dacs_strerror(err)
                  << "; Cell process (" << accel_filename << ") "
                  << "terminated with exit code " << exit_status
                  << std::endl;

    // Release the reserved Cell.
    err = dacs_release_de_list(1, &de_id);
    if (err != DACS_SUCCESS)
        std::cerr << dacs_strerror(err) << std::endl;

    // Shut down DACS.
    err = dacs_exit();
    if (err != DACS_SUCCESS)
        std::cerr << dacs_strerror(err) << std::endl;

} // DACS_Device::~DACS_Device

//---------------------------------------------------------------------------//
/*! \brief Record the full canonical path to the accel-side binary.
 *
 * init tries to interpret the specified filename in three ways: as an
 * absolute pathname, as a pathname relative to the current directory, and as
 * a pathname relative to a default location (default_ppe_bindir).
 * DACS_Device will launch exactly one accel-side binary, so calling init
 * multiple times with different filenames triggers an exception.
 */
void DACS_Device::init(const std::string & fname)
{
    // Identify the current working directory.
    char curr_path[MAXPATHLEN];
    Insist(getcwd(curr_path, MAXPATHLEN) != NULL,
           verbose_error("getcwd failed: " + std::string(strerror(errno))));

    const std::string cwd(curr_path);

    const std::string fullname1 = cwd + "/" + fname;
    const std::string fullname2 = default_ppe_bindir + "/" + fname;

    // Interpret fname as an absolute pathname...
    std::string fullname = canonical_fname(fname);
    if (fullname.empty())
    {
        // ... as a pathname relative to the current directory...
        fullname = canonical_fname(fullname1);
        if (fullname.empty())
        {
            // ... and as a pathname relative to default_ppe_bindir.
            fullname = canonical_fname(fullname2);
            Insist(!fullname.empty(),
                   verbose_error("Couldn't stat " + fname +
                                 " or " + fullname1 +
                                 " or " + fullname2 + ": " +
                                 std::string(strerror(errno))));
        }
    }

    if (accel_filename.empty())
        accel_filename = fullname;
    else
        Insist(fullname == accel_filename,
               "init was called twice with different filenames");
}

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
//                 end of DACS_Device.cc
//---------------------------------------------------------------------------//
