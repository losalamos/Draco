//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/DACS_Device_Interface.cc
 * \author Gabriel M. Rockefeller
 * \date   Fri Jun 24 12:58:18 2011
 * \brief  DACS_Device_Interface implementation.
 * \note   Copyright (C) 2011-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
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

// The definition in the source file does not have extern "C" on it, but
// because it is the same f(int) as in the declaration in the header, the
// language rules state that the extern "C" from the declaration applies to it
// as well. The definition for f is thus generated with C linkage, that is,
// with the unmangled name "f". (www.glenmccl.com/ansi028.htm).

// ExternCFuncPtrVoidVoid is defined in DACS_External_Process.hh

int dacs_device_init(const char * const filename, const int len, ExternCFuncPtrVoidVoid f)
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
// end of DACS_Device_Interface.cc
//---------------------------------------------------------------------------//
