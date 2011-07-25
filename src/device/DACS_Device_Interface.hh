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

#include <dacs.h>

extern "C"
{
    int dacs_device_init(const char * const filename, const int len,
                         void (*f)());

    int dacs_device_get_de_id(de_id_t * const de);

    int dacs_device_get_pid(dacs_process_id_t * const pid);
}

#endif // device_DACS_Device_Interface_hh

//---------------------------------------------------------------------------//
//              end of device/DACS_Device_Interface.hh
//---------------------------------------------------------------------------//
