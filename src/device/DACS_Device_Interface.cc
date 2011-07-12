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

using namespace rtt_device;

int DACS_DEVICE_INIT(const char * const filename, const int len)
{
    DACS_Device::init(std::string(filename, len));

    // Return success.
    return 0;
}

int DACS_DEVICE_GET_DE_ID(de_id_t * const de)
{
    *de = DACS_Device::instance().get_de_id();

    // Return success.
    return 0;
}

int DACS_DEVICE_GET_PID(dacs_process_id_t * const pid)
{
    *pid = DACS_Device::instance().get_pid();

    // Return success.
    return 0;
}

//---------------------------------------------------------------------------//
//                 end of DACS_Device_Interface.cc
//---------------------------------------------------------------------------//
