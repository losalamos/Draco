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
#include <vector>
#include <cmath>

#include "ds++/Assert.hh"
#include "ds++/Release.hh"
#include "c4/ParallelUnitTest.hh"

#include "DACS_Device_Interface.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

/*! \brief Test instance without calling init first.
 *
 * DACS_Device needs to know the name of the accel-side binary to launch.
 * Invoking DACS_DEVICE_GET_DE_ID or DACS_DEVICE_GET_PID without first calling
 * DACS_DEVICE_INIT(filename) should trigger an exception.
 */
void tstNoInit(UnitTest &ut)
{
    bool caught = false;
    de_id_t de_id;

    try
    {
        DACS_DEVICE_GET_DE_ID(&de_id);
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
 * dacs_de_start.  If the binary (specified via DACS_DEVICE_INIT) can't be
 * found, DACS_Device should throw an exception with a useful error message.
 */
void tstNoAccelBinary(UnitTest &ut)
{
    bool caught = false;
    string filename("no_such_binary");
    
    try
    {
        DACS_DEVICE_INIT(filename.c_str(), filename.size());
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
 * DACS_DEVICE_INIT must be called at least once.  Calling it multiple times
 * with different filenames is a sign of some kind of confusion, and
 * DACS_Device should throw an exception.  Calling it multiple times with the
 * same filename should be allowed, though.
 */
void tstDoubleInit(UnitTest &ut)
{
    bool caught = false;
    string f1("dacs_device_ppe_exe");
    string f2("dacs_device_ppe2_exe");

    try
    {
        // Call DACS_DEVICE_INIT twice with different filenames.
        DACS_DEVICE_INIT(f1.c_str(), f1.size());
        DACS_DEVICE_INIT(f2.c_str(), f2.size());
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

    // Call DACS_DEVICE_INIT again with the first filename.  This should not
    // trigger an exception.
    DACS_DEVICE_INIT(f1.c_str(), f1.size());
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
    string filename("dacs_device_ppe_exe");

    DACS_DEVICE_INIT(filename.c_str(), filename.size());

    de_id_t de_id;
    DACS_DEVICE_GET_DE_ID(&de_id);
    if (de_id == 0) ut.failure("de_id == 0");
    else ut.passes("de_id != 0");

    dacs_process_id_t pid;
    DACS_DEVICE_GET_PID(&pid);
    if (pid == 0) ut.failure("pid == 0");
    else ut.passes("pid != 0");

    de_id_t de_id2;
    DACS_DEVICE_GET_DE_ID(&de_id2);
    if (de_id2 != de_id)
        ut.failure("de_id changed");
    else
        ut.passes("de_id hasn't changed");

    dacs_process_id_t pid2;
    DACS_DEVICE_GET_PID(&pid2);
    if (pid2 != pid)
        ut.failure("pid changed");
    else
        ut.passes("pid hasn't changed");
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
