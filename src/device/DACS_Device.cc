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

#include "ds++/Assert.hh"

#include <cstring>
#include <iostream>
#include <sstream>
#include <stdint.h>
#include <sys/param.h>

namespace rtt_device
{

DACS_Device::SP_Proc DACS_Device::process;

static std::string verbose_error(std::string const & message);

//---------------------------------------------------------------------------//
/*! \brief DACS_Device default constructor.
 *
 * DACS_Device is a singleton, so the constructor is called at most once, the
 * first time any other code invokes DACS_Device::instance.  The constructor
 * reserves a Cell (via dacs_reserve_children) and launches the accel-side
 * binary by launching the DACS_Process that was passed to init.  If no
 * accel-side process has been specified (via DACS_Device::init), the
 * constructor throws an exception.
 */
DACS_Device::DACS_Device() : de_id(0)
{
    Insist(process, "No accel process specified");

    // Get the number of available Cell BEs.
    uint32_t num_children(0);
    DACS_ERR_T err = dacs_get_num_avail_children(DACS_DE_CBE, &num_children);
    Insist(err == DACS_SUCCESS, verbose_error(dacs_strerror(err)));

    Insist(num_children > 0, "No child CBE available");

    // Reserve one Cell BE.
    uint32_t num_reserve(1);
    err = dacs_reserve_children(DACS_DE_CBE, &num_reserve, &de_id);
    Insist(err == DACS_SUCCESS, verbose_error(dacs_strerror(err)));

    Insist(num_reserve == 1, "Failed to reserve children");

    // Start the accel-side process.
    process->start(de_id);

} // DACS_Device::DACS_Device

//---------------------------------------------------------------------------//
/*! \brief DACS_Device destructor.
 *
 * DACS_Device is a singleton; the destructor is called at program
 * termination.  The destructor tells the DACS_Process to shut down, waits for
 * the accel-side process to complete, releases the reserved Cell, and then
 * calls dacs_exit.  Throwing an exception inside the destructor would
 * interfere with stack unwinding already in progress (from, for example, an
 * earlier exception), so errors are reported to stderr.
 */
DACS_Device::~DACS_Device()
{
    // Tell the accel-side process to stop.
    process->stop();

    // Wait for the accel-side process to exit.
    int32_t exit_status;
    DACS_ERR_T err = dacs_de_wait(de_id,
                                  process->get_pid(),
                                  &exit_status);
    if (err != DACS_STS_PROC_FINISHED)
        std::cerr << dacs_strerror(err)
                  << "; Cell process (" << process->get_filename() << ") "
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
/*! \brief Initialize DACS and associate with a DACS_Process.
 *
 * init saves the specified DACS_Process in a static member variable, so that
 * the DACS_Device constructor can launch the accel-side process, and so that
 * the singleton can ask the process for its pid and tell it to stop during
 * host-side program termination.  init also initializes DACS, the first time
 * it's called.  DACS_Device will launch exactly one accel-side binary, so
 * calling init multiple times with different DACS_Processes that point to
 * different filenames triggers an exception.
 */
void DACS_Device::init(SP_Proc proc)
{
    // Initialize DACS.
    static bool dacs_is_initialized(false);

    if (!dacs_is_initialized)
    {
        DACS_ERR_T err = dacs_init(DACS_INIT_FLAGS_NONE);
        Insist(err == DACS_SUCCESS, verbose_error(dacs_strerror(err)));
        dacs_is_initialized = true;
    }

    if (!process)
        process = proc;
    else
        Insist(*process == *proc,
               "init was called twice with different filenames");
} // DACS_Device::init

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
