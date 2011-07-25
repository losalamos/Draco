//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/tstDACS_Device_Process.cc.cc
 * \author Gabriel M. Rockefeller
 * \date   Mon Jul 25 13:13:38 2011
 * \brief  DACS_Device + DACS_Process tests.
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>

#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "c4/ParallelUnitTest.hh"

#include "DACS_Device_Interface.hh"
#include "DACS_External_Process.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_device;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

extern "C"
{
    void dacs_process_shutdown()
    {
        de_id_t de_id;
        dacs_process_id_t pid;

        uint32_t command = 30306174;
        dacs_wid_t wid;

        // Reserve a wid.
        DACS_ERR_T err = dacs_wid_reserve(&wid);
        if (err != DACS_SUCCESS)
            std::cerr << "ERROR: " << dacs_strerror(err) << std::endl;

        // Get the de_id of the device, and the pid of the accel-side process.
        dacs_device_get_de_id(&de_id);
        dacs_device_get_pid(&pid);

        // Send the "stop" command expected by dacs_wait_for_cmd_ppe_exe.
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

        // Wait for the send to complete.
        err = dacs_wait(wid);
        if (err != DACS_WID_READY)
            std::cerr << "ERROR: " << verbose_error(dacs_strerror(err))
                      << std::endl;

        // Release the wid.
        err = dacs_wid_release(&wid);
        if (err != DACS_SUCCESS)
            std::cerr << "ERROR: " << verbose_error(dacs_strerror(err))
                      << std::endl;
    }
}

//---------------------------------------------------------------------------//
/*! \brief Exercise the DACS_Device destructor.
 *
 * This test exercises the accel-side stop-message mechanism through
 * DACS_External_Process.  dacs_device_init accepts a pointer to a function
 * that should be called to stop the accel-side process.  If the DACS_Device
 * destructor doesn't tell the accel-side process to stop, this test will
 * hang.
 */
void tstStartStop(UnitTest &ut)
{
    string filename("dacs_wait_for_cmd_ppe_exe");
    
    int rc = dacs_device_init(filename.c_str(), filename.length(),
                              dacs_process_shutdown);
    if (rc != 0) ut.failure("rc != 0");
    else ut.passes("rc == 0");

    de_id_t de_id;
    rc = dacs_device_get_de_id(&de_id);
    if (de_id == 0) ut.failure("de_id == 0");
    else ut.passes("de_id != 0");
    if (rc != 0) ut.failure("rc != 0");
    else ut.passes("rc == 0");

    dacs_process_id_t pid;
    rc = dacs_device_get_pid(&pid);
    if (pid == 0) ut.failure("pid == 0");
    else ut.passes("pid != 0");
    if (rc != 0) ut.failure("rc != 0");
    else ut.passes("rc == 0");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, release);
    try
    {
        tstStartStop(ut);
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing tstDACS_Device_Process.cc, " 
             << err.what()
             << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        cout << "ERROR: While testing tstDACS_Device_Process.cc, " 
             << "An unknown exception was thrown."
             << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstDACS_Device_Process.cc.cc
//---------------------------------------------------------------------------//
