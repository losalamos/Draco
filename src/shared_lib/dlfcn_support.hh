//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   shared_lib/dlfcn_supported.hh
 * \author Rob Lowrie
 * \date   Thu Apr 15 20:44:39 2004
 * \brief  Logic for dlfcn.h system header.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved.
 */
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

#ifndef rtt_shared_lib_dlfcn_supported_hh
#define rtt_shared_lib_dlfcn_supported_hh

// This header isolates the use of the macro NO_DLOPEN for Shared_Lib.
// It defines Shared_Lib::is_supported(), and whether the platform is
// supported or not, defines the system functions.

#ifdef USE_DLOPEN

// ... then dynamic loading of shared libraries is supported.

// Load the system header.

#include <dlfcn.h>

bool rtt_shared_lib::Shared_Lib::is_supported() { return true; }

#else

// ... then dynamic loading of shared libraries is unsupported.

bool rtt_shared_lib::Shared_Lib::is_supported() { return false; }

// Define dlfcn function stubs and constants so that Shared_Lib will compile.

static const int RTLD_LAZY = 1; // or whatever.

void *dlopen(const char * /*filename*/, int /*flag*/) {
  Insist(0, "Serious Shared_Lib error.");
  void *dummy(0);
  return dummy;
}

const char *dlerror(void) {
  Insist(0, "Serious Shared_Lib error.");
  return "dummy";
}

void *dlsym(void * /*handle*/, const char * /*symbol*/) {
  Insist(0, "Serious Shared_Lib error.");
  void *dummy(0);
  return dummy;
}

int dlclose(void * /*handle*/) {
  Insist(0, "Serious Shared_Lib error.");
  return 1;
}

#endif // USE_DLOPEN

#endif // rtt_shared_lib_dlfcn_supported_hh

//----------------------------------------------------------------------------//
// end of shared_lib/dlfcn_supported.hh
//----------------------------------------------------------------------------//
