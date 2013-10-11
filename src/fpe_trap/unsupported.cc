//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/unsupported.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 16:52:05 2005
 * \brief  Default implementation of fpe_trap functions.
 * \note   Copyright (C) 2005-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "fpe_trap/config.h"
#include "fpe_trap.hh"

#ifdef FPETRAP_UNSUPPORTED

DLL_PUBLIC bool rtt_fpe_trap::enable_fpe( bool /*AbortWithInsist*/ )
{
    return false;
}

#endif // FPETRAP_UNSUPPORTED

//---------------------------------------------------------------------------//
// end of unsupported.cc
//---------------------------------------------------------------------------//
