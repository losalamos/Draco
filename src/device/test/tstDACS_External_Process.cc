//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/tstDACS_External_Process.cc
 * \author Gabriel M. Rockefeller
 * \date   Thu Jul 14 16:39:58 2011
 * \brief  DACS_External_Process tests.
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>

#include <dacs.h>
#include <stdint.h>

#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "c4/ParallelUnitTest.hh"

#include "DACS_External_Process.hh"

#include "device/config.h"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_device;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

static de_id_t de_id = 0;
static dacs_process_id_t pid = 0;

extern "C"
{
    void dacs_process_shutdown()
    {
        uint32_t command = 30306174;
        dacs_wid_t wid;

        DACS_ERR_T err = dacs_wid_reserve(&wid);
        if (err != DACS_SUCCESS)
            std::cerr << "ERROR: " << dacs_strerror(err) << std::endl;

        err = dacs_send(&command, 
                        sizeof(uint32_t), 
                        de_id,
                        pid,
                        1, 
                        wid,
                        DACS_BYTE_SWAP_DISABLE);
        if (err != DACS_SUCCESS)
            std::cerr << "ERROR: " << verbose_error(dacs_strerror(err))
                      << std::endl;

        err = dacs_wait(wid);
        if (err != DACS_WID_READY)
            std::cerr << "ERROR: " << verbose_error(dacs_strerror(err))
                      << std::endl;

        err = dacs_wid_release(&wid);
        if (err != DACS_SUCCESS)
            std::cerr << "ERROR: " << verbose_error(dacs_strerror(err))
                      << std::endl;
    }
}

//---------------------------------------------------------------------------//
/*! \brief Test DACS_Process when the accel-side binary doesn't exist.
 *
 * The DACS_Process constructor checks for the existence of the accel-side
 * binary and should throw an exception with a useful error message if it
 * doesn't exist.
 */
void tstNoSuchBinary(UnitTest &ut)
{
    bool caught(false);

    try
    {
        DACS_External_Process d("no_such_binary", NULL);
    }
    catch (assertion &err)
    {
        caught = true;
        ostringstream message;
        message << "Good, caught the following assertion, " << err.what();
        ut.passes(message.str());
    }
    if (!caught)
    {
        ut.failure("Failed to catch expected assertion");
    }
}

//---------------------------------------------------------------------------//
/*! \brief Test the get_filename accessor.
 *
 * DACS_Process holds onto the full canonical pathname to the accel-side
 * binary.  Two DACS_Process instances pointing to the same accel-side binary
 * should compare as equal.
 */
void tstGetFilename(UnitTest &ut)
{
    const std::string filename("dacs_noop_ppe_exe");

    DACS_External_Process d(test_ppe_bindir + "/" + filename, NULL);
    std::string f = d.get_filename();

    if (f == test_ppe_bindir + "/" + filename)
        ut.passes("Filenames match");
    else
        ut.failure("Filenames don't match");

    const std::string fullpath(test_ppe_bindir + "/" + filename);

    DACS_External_Process d2(fullpath, NULL);
    f = d2.get_filename();

    if (f == test_ppe_bindir + "/" + filename)
        ut.passes("Filenames match");
    else
        ut.failure("Filenames don't match");

    if (d != d2) ut.failure("DACS_Processes aren't equal");
    else if (d == d2) ut.passes("DACS_Processes are equal");
    else ut.failure("Neither equal nor unequal?");
}

//---------------------------------------------------------------------------//
/*! \brief Test the start and stop methods.
 *
 * The start method should launch the accel-side binary.  get_pid should
 * return the pid of the accel-side process.  The stop method should invoke
 * the function specified during construction of DACS_External_Process to
 * signal the accel-side process to terminate; if the signal is not sent,
 * dacs_wait_for_cmd_ppe_exe won't exit, dacs_de_wait will block forever, and
 * this test will appear to hang.
 */
void tstStartStop(UnitTest &ut)
{
    // Initialize DACS.
    DACS_ERR_T err = dacs_init(DACS_INIT_FLAGS_NONE);

    if (err != DACS_SUCCESS)
        ut.failure(verbose_error(dacs_strerror(err)));

    // Get the number of available Cell BEs.
    uint32_t num_children(0);
    err = dacs_get_num_avail_children(DACS_DE_CBE, &num_children);

    if (err != DACS_SUCCESS)
        ut.failure(verbose_error(dacs_strerror(err)));

    if (num_children == 0)
        ut.failure("No child CBE available");

    // Reserve one Cell BE.
    uint32_t num_reserve(1);
    err = dacs_reserve_children(DACS_DE_CBE, &num_reserve, &de_id);

    if (err != DACS_SUCCESS)
        ut.failure(verbose_error(dacs_strerror(err)));

    if (num_reserve != 1)
        ut.failure("Failed to reserve children");

    // Start the accel-side process.
    DACS_External_Process d(test_ppe_bindir + "/dacs_wait_for_cmd_ppe_exe",
                            dacs_process_shutdown);

    d.start(de_id);

    // Create a group containing the host and the accel-side process.
    dacs_group_t group;
    err = dacs_group_init(&group, 0u);
    Insist(err == DACS_SUCCESS, verbose_error(dacs_strerror(err)));

    err = dacs_group_add_member(de_id, d.get_pid(), group);
    Insist(err == DACS_SUCCESS, verbose_error(dacs_strerror(err)));

    err = dacs_group_add_member(DACS_DE_SELF, DACS_PID_SELF, group);
    Insist(err == DACS_SUCCESS, verbose_error(dacs_strerror(err)));

    err = dacs_group_close(group);
    Insist(err == DACS_SUCCESS, verbose_error(dacs_strerror(err)));

    pid = d.get_pid();
    if (pid == 0) ut.failure("pid == 0");
    else ut.passes ("pid != 0");

    // Tell the accel-side process to stop.
    d.stop();

    // Destroy the host-accel group.
    err = dacs_group_destroy(&group);
    Insist(err == DACS_SUCCESS, verbose_error(dacs_strerror(err)));

    // Wait for the accel-side process to exit.
    int32_t exit_status;
    err = dacs_de_wait(de_id, pid, &exit_status);

    if (err != DACS_STS_PROC_FINISHED)
        ut.failure(verbose_error(dacs_strerror(err)));
           
    // Release the reserved Cell.
    err = dacs_release_de_list(1, &de_id);

    if (err != DACS_SUCCESS)
        ut.failure(verbose_error(dacs_strerror(err)));

    // Shut down DACS.
    err = dacs_exit();

    if (err != DACS_SUCCESS)
        ut.failure(verbose_error(dacs_strerror(err)));

    ut.passes("Cell process stopped successfully");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, release);
    try
    {
        tstNoSuchBinary(ut);
        tstGetFilename(ut);
        tstStartStop(ut);
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing tstDACS_External_Process, " 
             << err.what()
             << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        cout << "ERROR: While testing tstDACS_External_Process, " 
             << "An unknown exception was thrown."
             << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstDACS_External_Process.cc
//---------------------------------------------------------------------------//
