//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/darwin_ppc.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 16:52:05 2005
 * \brief  Darwin/PPC implementation of fpe_trap functions.
 *
 * Copyright (C) 2004-2013  Los Alamos National Security, LLC.
 *               All rights reserved.
 * Copyright (C) 1994-2001  K. Scott Hunziker.
 * Copyright (C) 1990-1994  The Boeing Company.
 *
 * See COPYING file for more copyright information.  This code is based
 * substantially on fpe/darwin.c from algae-4.3.6, which is available at
 * http://algae.sourceforge.net/.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <fpe_trap/config.h>

#ifdef FPETRAP_DARWIN_PPC

#include <iostream>
#include <ds++/Assert.hh>
#include <mach/mach.h>

// Local functions

namespace
{

/*
 * On Mach, we need a mindbogglingly complex setup for floating point errors.
 * Not the least of the hassles is that we have to do the whole thing from
 * a different thread.
 */
void* fpe_enabler(void *parent)
{
    mach_port_t		   victim = (mach_port_t)parent;
    mach_msg_type_number_t count;

    ppc_thread_state_t	ts;
    ppc_float_state_t	fs;

    /* First enable the right FP exception conditions */
    count = PPC_FLOAT_STATE_COUNT;
    thread_get_state(victim, PPC_FLOAT_STATE, (thread_state_t)&fs, &count);
    /* Enable VE OE ZE, Disable UE XE */
    fs.fpscr = (fs.fpscr & ~0x1FFFFF28) | 0x0D0;
    thread_set_state(victim, PPC_FLOAT_STATE, (thread_state_t)&fs, count);
    
    /* Now, enable FP exceptions as such */
    count = PPC_THREAD_STATE_COUNT;
    thread_get_state(victim, PPC_THREAD_STATE, (thread_state_t)&ts, &count);
    /* Set FE0 = FE1 = 1 */
    ts.srr1 |= 0x0900;
    thread_set_state(victim, PPC_THREAD_STATE, (thread_state_t)&ts, count);
    
    return 0;
}

static void catch_sigfpe(int sig)
{
    Insist(0, "Floating point exception caught by fpe_trap.")
}

} // end of namespace

//---------------------------------------------------------------------------------------//
namespace rtt_fpe_trap
{
//---------------------------------------------------------------------------------------//
//!  Enable trapping of floating point errors.
bool fpe_trap::enable(void)
{
    pthread_t enabler;
    void *mts = reinterpret_cast<void *>(mach_thread_self());
    pthread_create(&enabler, NULL, fpe_enabler, mts);
    pthread_join(enabler, NULL);

    if( this->abortWithInsist)
        signal(SIGFPE, catch_sigfpe);

     // Toggle the state.
    fpeTrappingActive = true;
    return fpeTrappingActive;
}

//---------------------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable(void)
{
    // (void)feenableexcept( 0x00 );
    Insist(0,"Please update darwin_ppc.cc to provide instructions for disabling fpe traps.");
    // fpeTrappingActive=false;
    return;
}

} // end namespace rtt_shared_lib

#endif // FPETRAP_DARWIN_PPC

//---------------------------------------------------------------------------//
// end of darwin_ppc.cc
//---------------------------------------------------------------------------//
