//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/linux_x86.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 16:52:05 2005
 * \brief  Linux/X86 implementation of fpe_trap functions.
 *
 * Copyright (C) 2005-2013 Los Alamos National Security, LLC.
 *               All rights reserved.
 * Copyright (C) 1994-2001  K. Scott Hunziker.
 * Copyright (C) 1990-1994  The Boeing Company.
 *
 * See COPYING file for more copyright information.  This code is based
 * substantially on fpe/i686-pc-linux-gnu.c from algae-4.3.6, which is
 * available at http://algae.sourceforge.net/.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "fpe_trap/config.h"
#include "fpe_trap.hh"

#ifdef FPETRAP_LINUX_X86

#include <ds++/Assert.hh>
#include <iostream>
#include <string>
#include <signal.h>
#include <fenv.h>

/* Signal handler for floating point exceptions. */
extern "C"
{

static void
catch_sigfpe (int sig, siginfo_t *code, void * /*v*/)
{
    std::cout << "(fpe_trap/linux_x86.cc) A SIGFPE was detected!"
              << std::endl;
        
    std::string mesg;
    if (sig != SIGFPE)
    {
        mesg = "Floating point exception problem.";
    }
    else
    {
        switch (code->si_code)
        {
            case FPE_INTDIV:
                mesg = "Integer divide by zero.";
                break;
            case FPE_INTOVF:
                mesg = "Integer overflow.";
                break;
            case FPE_FLTDIV:
                mesg = "Floating point divide by zero.";
                break;
            case FPE_FLTOVF:
                mesg = "Floating point overflow.";
                break;
            case FPE_FLTUND:
                mesg = "Floating point underflow.";
                break;
            case FPE_FLTRES:
                mesg = "Floating point inexact result.";
                break;
            case FPE_FLTINV:
                mesg = "Invalid floating point operation.";
                break;
            case FPE_FLTSUB:
                mesg = "Floating point subscript out of range.";
                break;
            default:
                mesg = "Unknown floating point exception.";
                break;
        }
    }
    
    Insist(0, mesg);
}

} // end of extern "C"

namespace rtt_fpe_trap
{

//---------------------------------------------------------------------------------------//
//!  Enable trapping of floating point errors.
bool fpe_trap::enable(void)
{
    struct sigaction act;

    // Choose to use Draco's DbC Insist.  If set to false, the compiler should
    // print a stack trace instead of the pretty print message defined above
    // in catch_sigfpe.
    if( this->abortWithInsist )
        act.sa_sigaction = catch_sigfpe; /* the signal handler       */

    sigemptyset(&(act.sa_mask));         /* no other signals blocked */
    act.sa_flags = SA_SIGINFO;           /* want 3 args for handler  */

    // specify handler
    Insist(! sigaction(SIGFPE, &act, NULL),
           "Unable to set floating point handler.");

    // The feenableexcept function is new for glibc 2.2.  See its description
    // in the man page for fenv(3). 

    (void)feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

    // Toggle the state.
    fpeTrappingActive = true;
                        
    // ... if you want to catch "everything" :
    //(void)feenableexcept(FE_ALL_EXCEPT);
    
    // Also available (taken from fenv.h and bits/fenv.h):
    // FE_INEXACT               inexact result
    // FE_DIVBYZERO             division by zero
    // FE_UNDERFLOW             result not representable due to underflow
    // FE_OVERFLOW              result not representable due to overflow
    // FE_INVALID               invalid operation

    // FE_ALL_EXCEPT            bitwise OR of all supported exceptions

    // The next macros are defined iff the appropriate rounding mode is
    // supported by the implementation.
    // FE_TONEAREST             round to nearest
    // FE_UPWARD                round toward +Inf
    // FE_DOWNWARD              round toward -Inf
    // FE_TOWARDZERO            round toward 0
    
    return fpeTrappingActive;
}

//---------------------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable(void)
{
    (void)feenableexcept( 0x00 );
    fpeTrappingActive=false;
    return;
}

} // end namespace rtt_shared_lib

#endif // FPETRAP_LINUX_X86

//---------------------------------------------------------------------------//
// end of linux_x86.cc
//---------------------------------------------------------------------------//
