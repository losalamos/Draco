//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/linux_x86.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 16:52:05 2005
 * \brief  Linux/X86 implementation of fpe_trap functions.
 *
 * Copyright 2004 The Regents of the University of California.
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

#include <fpe_trap/config.h>

#ifdef FPETRAP_LINUX_X86

#include <iostream>
#include <string>
#include <ds++/Assert.hh>

#include <signal.h>
#include <fenv.h>

// Local functions

namespace
{

/* Signal handler for floating point exceptions. */

static void
catch_sigfpe (int sig, siginfo_t *code, void * /*v*/)
{
//    using namespace std;
    
    std::string mesg;

//    cout << "In catch_sigfpe" << endl;
    
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

} // end of namespace

namespace rtt_fpe_trap
{

bool enable_fpe()
{
    struct sigaction act;

    act.sa_sigaction = catch_sigfpe;	/* the signal handler */
    sigemptyset(&(act.sa_mask));	/* no other signals blocked */
    act.sa_flags = SA_SIGINFO;		/* want 3 args for handler */

    /* specify handler */
    Insist(! sigaction(SIGFPE, &act, NULL),
           "Unable to set floating point handler.");

    /*
     * The feenableexcept function is new for glibc 2.2.  See its
     * description in the man page for fenv(3).
     */

    (void)feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    
    // ... if you want to catch "everything" :
    //(void)feenableexcept(FE_ALL_EXCEPT);
    // Also available:
    // FE_INEXACT
    // FE_UNDERFLOW

    return true;
}

} // end namespace rtt_shared_lib

#endif // FPETRAP_LINUX_X86

//---------------------------------------------------------------------------//
//                 end of linux_x86.cc
//---------------------------------------------------------------------------//
