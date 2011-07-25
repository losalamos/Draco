//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/DACS_Device_Interface.cc
 * \author Gabriel M. Rockefeller
 * \date   Fri Jun 24 12:58:18 2011
 * \brief  DACS_Device_Interface implementation.
 * \note   Copyright (C) 2006 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DACS_Device_Interface.hh"
#include "DACS_Device.hh"
#include "DACS_External_Process.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_device;

int dacs_device_init(const char * const filename, const int len, void (*f)())
{
    SP<DACS_Process> p(new DACS_External_Process(string(filename, len), f));

    DACS_Device::init(p);

    // Return success.
    return 0;
}

int dacs_device_get_de_id(de_id_t * const de)
{
    *de = DACS_Device::instance().get_de_id();

    // Return success.
    return 0;
}

int dacs_device_get_pid(dacs_process_id_t * const pid)
{
    *pid = DACS_Device::instance().get_pid();

    // Return success.
    return 0;
}

//---------------------------------------------------------------------------//
//                 end of DACS_Device_Interface.cc
//---------------------------------------------------------------------------//
