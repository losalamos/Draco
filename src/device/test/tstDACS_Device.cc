//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/tstDACS_Device.cc
 * \author Gabriel M. Rockefeller
 * \date   Mon Jun 13 16:54:19 2011
 * \brief  DACS_Device tests.
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC
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

#include "DACS_Device.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_device;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

/*! \brief Test instance without calling init first.
 *
 * DACS_Device needs to know the name of the accel-side binary to launch.
 * Invoking instance without first calling init(filename) should trigger an
 * exception.
 */
void tstNoInit(UnitTest &ut)
{
    bool caught = false;

    try
    {
        DACS_Device::instance();
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
 * init must be called at least once.  Calling it multiple times with
 * different filenames is a sign of some kind of confusion, and init should
 * throw an exception.  Calling it multiple times with the same filename
 * should be allowed, though.
 */
void tstDoubleInit(UnitTest &ut)
{
    bool caught = false;

    try
    {
        // Call init twice with different filenames.
        DACS_Device::init("dacs_device_ppe_exe");
        DACS_Device::init("some_other_file");
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

    // Call init again with the first filename.  This should not trigger an
    // exception.
    DACS_Device::init("dacs_device_ppe_exe");
    ut.passes("Called init twice with the same filename");
}

//---------------------------------------------------------------------------//
/*! \brief Test instance and the de_id and pid accessors.
 *
 * instance will create the DACS_Device on first invocation.  The de_id and
 * child pid shouldn't change after the DACS_Device has been created.
 */
void tstDevice(UnitTest &ut)
{
    DACS_Device::init("dacs_device_ppe_exe");

    de_id_t de_id = DACS_Device::instance().get_de_id();
    if (de_id == 0) ut.failure("de_id == 0");
    else ut.passes("de_id != 0");

    dacs_process_id_t pid = DACS_Device::instance().get_pid();
    if (pid == 0) ut.failure("pid == 0");
    else ut.passes("pid != 0");

    if (DACS_Device::instance().get_de_id() != de_id)
        ut.failure("de_id changed");
    else
        ut.passes("de_id hasn't changed");

    if (DACS_Device::instance().get_pid() != pid)
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
        tstDoubleInit(ut);
        tstDevice(ut);
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing tstDACS_Device, " 
             << err.what()
             << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        cout << "ERROR: While testing tstDACS_Device, " 
             << "An unknown exception was thrown."
             << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstDACS_Device.cc
//---------------------------------------------------------------------------//
