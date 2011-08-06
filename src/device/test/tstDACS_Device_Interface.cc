//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/tstDACS_Device_Interface.cc
 * \author Gabriel M. Rockefeller
 * \date   Fri Jun 24 13:49:31 2011
 * \brief  DACS_Device_Interface tests.
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <sstream>

#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "c4/ParallelUnitTest.hh"

#include "device/config.h"

#include "DACS_Device_Interface.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_device;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

/*! \brief Test instance without calling init first.
 *
 * DACS_Device needs to know the name of the accel-side binary to launch.
 * Invoking dacs_device_get_de_id or dacs_device_get_pid without first calling
 * dacs_device_init(filename) should trigger an exception.
 */
void tstNoInit(UnitTest &ut)
{
    bool caught = false;
    de_id_t de_id;

    try
    {
        dacs_device_get_de_id(&de_id);
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
/*! \brief Test DACS_Device when the accel-side binary doesn't exist.
 *
 * DACS_Device needs the path to an accel-side binary to successfully call
 * dacs_de_start.  If the binary (specified via dacs_device_init) can't be
 * found, DACS_Device should throw an exception with a useful error message.
 */
void tstNoAccelBinary(UnitTest &ut)
{
    bool caught = false;
    string filename(test_ppe_bindir + "/no_such_binary");
    
    try
    {
        dacs_device_init(filename.c_str(), filename.size(), NULL);
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
/*! \brief Test multiple calls to init, with different or identical filenames.
 *
 * dacs_device_init must be called at least once.  Calling it multiple times
 * with different filenames is a sign of some kind of confusion, and
 * DACS_Device should throw an exception.  Calling it multiple times with the
 * same filename should be allowed, though.
 */
void tstDoubleInit(UnitTest &ut)
{
    bool caught = false;
    string f1(test_ppe_bindir + "/dacs_noop_ppe_exe");
    string f2(test_ppe_bindir + "/dacs_wait_for_cmd_ppe_exe");

    try
    {
        // Call dacs_device_init twice with different filenames.
        dacs_device_init(f1.c_str(), f1.size(), NULL);
        dacs_device_init(f2.c_str(), f2.size(), NULL);
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

    // Call dacs_device_init again with the first filename.  This should not
    // trigger an exception.
    dacs_device_init(f1.c_str(), f1.size(), NULL);
    ut.passes("Called init twice with the same filename");
}

//---------------------------------------------------------------------------//
/*! \brief Test the de_id and pid accessors.
 *
 * The accessors invoke DACS_Device::instance, which will create the
 * DACS_Device on first use.  The de_id and child pid shouldn't change after
 * the DACS_Device has been created.
 */
void tstDevice(UnitTest &ut)
{
    string filename(test_ppe_bindir + "/dacs_noop_ppe_exe");
    int rc;

    rc = dacs_device_init(filename.c_str(), filename.size(), NULL);
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

    de_id_t de_id2;
    rc = dacs_device_get_de_id(&de_id2);
    if (de_id2 != de_id)
        ut.failure("de_id changed");
    else
        ut.passes("de_id hasn't changed");
    if (rc != 0) ut.failure("rc != 0");
    else ut.passes("rc == 0");

    dacs_process_id_t pid2;
    rc = dacs_device_get_pid(&pid2);
    if (pid2 != pid)
        ut.failure("pid changed");
    else
        ut.passes("pid hasn't changed");
    if (rc != 0) ut.failure("rc != 0");
    else ut.passes("rc == 0");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, release);
    try
    {
        tstNoInit(ut);
        tstNoAccelBinary(ut);
        tstDoubleInit(ut);
        tstDevice(ut);
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing tstDACS_Device_Interface, " 
             << err.what()
             << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        cout << "ERROR: While testing tstDACS_Device_Interface, " 
             << "An unknown exception was thrown."
             << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstDACS_Device_Interface.cc
//---------------------------------------------------------------------------//
