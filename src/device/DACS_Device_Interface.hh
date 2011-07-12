//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   device/DACS_Device_Interface.hh
 * \author Gabriel M. Rockefeller
 * \date   Fri Jun 24 12:58:18 2011
 * \brief  DACS_Device_Interface, for use by non-C++ clients.
 * \note   Copyright (C) 2011 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef device_DACS_Device_Interface_hh
#define device_DACS_Device_Interface_hh

#include "device/config.h"

#include <dacs.h>

extern "C"
{
    int DACS_DEVICE_INIT(const char * const filename, const int len);

    int DACS_DEVICE_GET_DE_ID(de_id_t * const de);

    int DACS_DEVICE_GET_PID(dacs_process_id_t * const pid);
}

#endif // device_DACS_Device_Interface_hh

//---------------------------------------------------------------------------//
//              end of device/DACS_Device_Interface.hh
//---------------------------------------------------------------------------//
