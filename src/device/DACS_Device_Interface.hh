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
#include "DACS_External_Process.hh"

// Note: That these 3 functions will be generated with C linkage (no name
// mangling).  However, the function pointer will still be generated with C++
// linkage in the .cc file and with C linkage in this file (according to the
// ANSI standard). This causes link errors with some compilers (PGI) because
// the function signatures in the .cc and .hh files do not match.  The
// solution is to create and use a typdef for the extern function pointer
// signature.

extern "C"
{
    // As defined in DACS_External_Process.hh:
    // typedef void (*ExternCFuncPtrVoidVoid)(void);
    
    int dacs_device_init(const char * const filename,
                         const int len,
                         ExternCFuncPtrVoidVoid);

    int dacs_device_get_de_id(de_id_t * const de);

    int dacs_device_get_pid(dacs_process_id_t * const pid);
}

#endif // device_DACS_Device_Interface_hh

//---------------------------------------------------------------------------//
//              end of device/DACS_Device_Interface.hh
//---------------------------------------------------------------------------//
