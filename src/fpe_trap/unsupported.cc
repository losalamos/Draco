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


namespace rtt_fpe_trap
{

//---------------------------------------------------------------------------------------//
//!  Enable trapping of floating point errors.
bool fpe_trap::enable(void)
{
    // Toggle the state.
    fpeTrappingActive = true;
    return fpeTrappingActive;
}
//---------------------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable(void)
{
    return;
}

#endif // FPETRAP_UNSUPPORTED

//---------------------------------------------------------------------------//
// end of unsupported.cc
//---------------------------------------------------------------------------//
