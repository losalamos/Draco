//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/osf_alpha.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 16:52:05 2005
 * \brief  OSF alpha implementation of fpe_trap functions.
 *
 * Copyright (C) 2004-2013  Los Alamos National Security, LLC.
 * Copyright (C) 1994-2001  K. Scott Hunziker.
 * Copyright (C) 1990-1994  The Boeing Company.
 *
 * See COPYING file for more copyright information.  This code is based
 * substantially on fpe/alpha-dec-osf3.0.c from algae-4.3.6, which is
 * available at http://algae.sourceforge.net/.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <fpe_trap/config.h>

#ifdef FPETRAP_OSF_ALPHA

#include <iostream>
#include <string>
#include <ds++/Assert.hh>

#include <signal.h>
#include <machine/fpu.h>

// Local functions

extern "C"
{

/* Signal handler for floating point exceptions. */

static void
catch_sigfpe(int sig)
{
    std::string mesg = "Floating point exception";

    // decipher sig later...

    Insist(0, mesg);
}

} // end of namespace

namespace rtt_fpe_trap
{

bool enable_fpe( bool abortWithInsist )
{
    unsigned long csr = ieee_get_fp_control();
    csr |= IEEE_TRAP_ENABLE_INV | IEEE_TRAP_ENABLE_DZE | IEEE_TRAP_ENABLE_OVF;
    ieee_set_fp_control(csr);

    if( abortWithInsist )
        signal(SIGFPE, catch_sigfpe);

    return true;
}

} // end namespace rtt_shared_lib

#endif // FPETRAP_OSF_ALPHA

//---------------------------------------------------------------------------//
// end of osf_alpha.cc
//---------------------------------------------------------------------------//
