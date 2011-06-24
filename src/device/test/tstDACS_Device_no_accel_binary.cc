//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/test/tstDACS_Device_no_accel_binary.cc
 * \author Gabriel M. Rockefeller
 * \date   Thu Jun 23 17:59:23 2011
 * \brief  DACS_Device tests, when the accel-side binary doesn't exist.
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

/*! \brief Test DACS_Device when the accel-side binary doesn't exist.
 *
 * DACS_Device needs the path to an accel-side binary to successfully call
 * dacs_de_start.  If the binary (specified via DACS_Device::init) can't be
 * found, DACS_Device should throw an exception with a useful error message.
 *
 * (Because DACS_Device doesn't allow the accel-side binary name to be "reset"
 * by subsequent calls to DACS_Device::init, this test needs to live in its
 * own binary.)
 */
void tstNoAccelBinary(UnitTest &ut)
{
    DACS_Device::init("no_such_binary");

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

int main(int argc, char *argv[])
{
    rtt_c4::ParallelUnitTest ut(argc, argv, release);
    try
    {
        tstNoAccelBinary(ut);
    }
    catch (exception &err)
    {
        cout << "ERROR: While testing tstDACS_Device_no_accel_binary, " 
             << err.what()
             << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        cout << "ERROR: While testing tstDACS_Device_no_accel_binary, " 
             << "An unknown exception was thrown."
             << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of tstDACS_Device_no_accel_binary.cc
//---------------------------------------------------------------------------//
