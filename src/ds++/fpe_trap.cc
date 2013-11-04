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

#include "fpe_trap.hh"
#include "Assert.hh"
#include <iostream>
#include <sstream>
#include <string>

//---------------------------------------------------------------------------//
// Linux_x86
//---------------------------------------------------------------------------//
#ifdef FPETRAP_LINUX_X86

#include <signal.h>
#include <fenv.h>

// Print a demangled stack backtrace of the caller function to FILE* out.
// http://blog.aplikacja.info/2010/12/backtraces-for-c/
#include <cxxabi.h>
#include <execinfo.h> // backtrace
#include <stdio.h>  //snprintf
#include <stdlib.h>
#include <string.h>
#include <ucontext.h>
#include <unistd.h>
// #include <Qt/qapplication.h>

/* Signal handler for floating point exceptions. */
extern "C"
{

// Print a demangled stack backtrace of the caller function to FILE* out.
// http://blog.aplikacja.info/2010/12/backtraces-for-c/
std::string
print_stacktrace( std::string const & error_name )
{
    // max size of stack backtrace.
    unsigned const max_frames(63);

    // store the message here.
    std::ostringstream msg;

    unsigned const linkname_size(512);
    char linkname[linkname_size]; 
    char buf[linkname_size];
    pid_t pid;
    int ret;

    /* Get our PID and build the name of the link in /proc */
    pid = getpid();
    snprintf(linkname, sizeof(linkname), "/proc/%i/exe", pid);
        
    /* Now read the symbolic link */
    ret = readlink(linkname, buf, linkname_size);
    buf[ret] = 0;

    msg << "\nStack trace:"
        << "\n  Signaling error: " << error_name
        << "\n  Process " << buf << "\n\n";
    
    // storage array for stack trace address data
    void * addrlist[max_frames+1];
    
    // retrieve current stack addresses
    int addrlen = backtrace(addrlist, sizeof(addrlist) / sizeof(void*));

    if (addrlen == 0)
    {
        msg << "  \n";
        return msg.str();
    }

    // resolve addresses into strings containing "filename(function+address)",
    // this array must be free()-ed
    char** symbollist = backtrace_symbols(addrlist, addrlen);

    // allocate string which will be filled with the demangled function name
    size_t funcnamesize = 256;
    char* funcname = (char*)malloc(funcnamesize);

    // iterate over the returned symbol lines. skip first two,
    // (addresses of this function and handler)
    for (int i = 2; i < addrlen; i++)
    {
	char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

	// find parentheses and +address offset surrounding the mangled name:
	// ./module(function+0x15c) [0x8048a6d]
	for (char *p = symbollist[i]; *p; ++p)
	{
	    if (*p == '(')
		begin_name = p;
	    else if (*p == '+')
		begin_offset = p;
	    else if (*p == ')' && begin_offset) {
		end_offset = p;
		break;
	    }
	}

	if (begin_name && begin_offset && end_offset
	    && begin_name < begin_offset)
	{
	    *begin_name++ = '\0';
	    *begin_offset++ = '\0';
	    *end_offset = '\0';

	    // mangled name is now in [begin_name, begin_offset) and caller
	    // offset in [begin_offset, end_offset). now apply
	    // __cxa_demangle():

	    int status;
	    char* ret = abi::__cxa_demangle(begin_name,
					    funcname, &funcnamesize, &status);
	    if (status == 0)
            {
		funcname = ret; // use possibly realloc()-ed string
                msg << "  (PID:" << pid << ") "
                    << symbollist[i] << " : "
                    << funcname << "()+" << begin_offset << "\n";
	    }
	    else
            {
		// demangling failed. Output function name as a C function with
		// no arguments.
                msg << "  (PID:" << pid << ") "
                    << symbollist[i] << " : "
                    << begin_name << "()+" << begin_offset << "\n";
	    }
	}
	else
	{
	    // couldn't parse the line? print the whole line.
            msg << "  (PID:" << pid << ") "
                << symbollist[i] << " : ??\n";
	}
    }

    free(funcname);
    free(symbollist);

    // fprintf(out, "stack trace END (PID:%d)\n", pid);
    msg << "Stack trace: END (PID:" << pid << ")\n";
    return msg.str();
}


//---------------------------------------------------------------------------//
static void
catch_sigfpe (int sig, siginfo_t *psSiginfo, void * /*psContext*/) 
{    
    // generate a message:
    std::string error_type;

    if (sig != SIGFPE)
    {
        error_type = "Floating point exception problem.";
    }
    else
    {
        switch (psSiginfo->si_code)
        {
            case FPE_INTDIV:
                error_type = "SIGFPE (Integer divide by zero)";
                break;
            case FPE_INTOVF:
                error_type = "SIGFPE (Integer overflow)";
                break;
            case FPE_FLTDIV:
                error_type = "SIGFPE (Floating point divide by zero)";
                break;
            case FPE_FLTOVF:
                error_type = "SIGFPE (Floating point overflow)";
                break;
            case FPE_FLTUND:
                error_type = "SIGFPE (Floating point underflow)";
                break;
            case FPE_FLTRES:
                error_type = "SIGFPE (Floating point inexact result)";
                break;
            case FPE_FLTINV:
                error_type = "SIGFPE (Invalid floating point operation)";
                break;
            case FPE_FLTSUB:
                error_type = "SIGFPE (Floating point subscript out of range)";
                break;
            default:
                error_type = "SIGFPE (Unknown floating point exception)";
                break;
        }
    }
    Insist(false, print_stacktrace(error_type));
}

} // end of extern "C"

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable(void)
{
    (void)fedisableexcept( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    fpeTrappingActive=false;
    return;
}

} // end namespace rtt_shared_lib

#endif // FPETRAP_LINUX_X86

//---------------------------------------------------------------------------//
// OSF_ALPHA
//---------------------------------------------------------------------------//
#ifdef FPETRAP_OSF_ALPHA

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

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
//!  Enable trapping of floating point errors.
bool fpe_trap::enable(void)
{
    unsigned long csr = ieee_get_fp_control();
    csr |= IEEE_TRAP_ENABLE_INV | IEEE_TRAP_ENABLE_DZE | IEEE_TRAP_ENABLE_OVF;
    ieee_set_fp_control(csr);

    if( abortWithInsist )
        signal(SIGFPE, catch_sigfpe);

    // Toggle the state.
    fpeTrappingActive = true;
    return fpeTrappingActive;
}

//---------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable(void)
{
    ieee_set_fp_control(0x00);
    fpeTrappingActive=false;
    return;
}

} // end namespace rtt_shared_lib

#endif // FPETRAP_OSF_ALPHA

//---------------------------------------------------------------------------//
// WINDOWS X86
//---------------------------------------------------------------------------//
#ifdef FPETRAP_WINDOWS_X86

#include <Windows.h>  // defines STATUS_FLOAT_...
#include <float.h>    // defines _controlfp_s

#pragma fenv_access (on)

/* Signal handler for floating point exceptions. */
extern "C" void trans_func( unsigned int u, PEXCEPTION_POINTERS pExp )
{
    std::cout << "(fpe_trap/windows_x86.cc) A SIGFPE was detected!"
              << std::endl;
    
    std::string mesg;
    switch (u)
    {
        case STATUS_INTEGER_DIVIDE_BY_ZERO:
            mesg = "Integer divide by zero.";
            break;
        case STATUS_INTEGER_OVERFLOW:
            mesg = "Integer overflow.";
            break;
        case STATUS_FLOAT_DIVIDE_BY_ZERO:
            mesg = "Floating point divide by zero.";
            break;
        case STATUS_FLOAT_OVERFLOW:
            mesg = "Floating point overflow.";
            break;
        case STATUS_FLOAT_UNDERFLOW:
            mesg = "Floating point underflow.";
            break;
        case STATUS_FLOAT_INEXACT_RESULT:
            mesg = "Floating point inexact result.";
            break;
        case STATUS_FLOAT_INVALID_OPERATION:
            mesg = "Invalid floating point operation.";
            break;
        default:
            mesg = "Unknown floating point exception.";
            break;
    }

    Insist(0, mesg);
}

namespace rtt_dsxx
{

// ----------------------------------------------------------------------------
// - http://stackoverflow.com/questions/2769814/how-do-i-use-try-catch-to-catch-floating-point-errors
// - See MSDN articles on fenv_access and _controlfp_s examples.
// ----------------------------------------------------------------------------
bool fpe_trap::enable( void )
{   
   // Allways call this before setting control words.
   _clearfp();

   // Read the current control words.
   unsigned int fp_control_word = _controlfp(0,0);

   // Set the exception masks off for exceptions that you want to trap.  When
   // a mask bit is set, the corresponding floating-point exception is 
   // blocked from being generated.

   fp_control_word &= ~( EM_INVALID | EM_ZERODIVIDE | EM_OVERFLOW );
   
   // Other options:
   // _EM_DENORMAL
   // _EM_UNDERFLOW
   // _EM_INEXACT
   
   // Update the control word with our changes
   // MCW_EM is Interrupt exception mask.
   _controlfp(fp_control_word,MCW_EM);

   if( this->abortWithInsist )
       _set_se_translator(trans_func);
       
   // Toggle the state.
    fpeTrappingActive = true;
   
   return fpeTrappingActive;
}

//---------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable(void)
{
   // Allways call this before setting control words.
   _clearfp();

   // Read the current control words.
   unsigned int fp_control_word = _controlfp(0,0);

   // Set the exception masks off for exceptions that you want to trap.  When
   // a mask bit is set, the corresponding floating-point exception is 
   // blocked from being generated.
   fp_control_word |= ( EM_INVALID | EM_ZERODIVIDE | EM_OVERFLOW );

   // Update the control word with our changes
   // MCW_EM is Interrupt exception mask.
   _controlfp(fp_control_word,MCW_EM);   
    
    fpeTrappingActive=false;
    return;
}

} // end namespace rtt_shared_lib

#endif // FPETRAP_WINDOWS_X86

//---------------------------------------------------------------------------//
// DARWIN INTEL
//---------------------------------------------------------------------------//
#ifdef FPETRAP_DARWIN_INTEL

#include <xmmintrin.h>

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
//! Enable trapping of floating point errors.
bool fpe_trap::enable(void)
{
    _mm_setcsr( _MM_MASK_MASK &~
                (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );    

    // This functionality is not currently implemented on the Mac.
    abortWithInsist = false;
    
    // Toggle the state.
    fpeTrappingActive = true;
    return fpeTrappingActive;
}
//---------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable(void)
{
    _mm_setcsr( _MM_MASK_MASK &~ 0x00 );  
    fpeTrappingActive=false;  
    return;
}

} // end namespace rtt_shared_lib

#endif // FPETRAP_DARWIN_INTEL

//---------------------------------------------------------------------------//
// DARWIN PPC
//---------------------------------------------------------------------------//
#ifdef FPETRAP_DARWIN_PPC

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

namespace rtt_dsxx
{
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

//---------------------------------------------------------------------------//
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
// UNSUPPORTED PLATFORMS
//---------------------------------------------------------------------------//
#ifdef FPETRAP_UNSUPPORTED

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
//!  Enable trapping of floating point errors.
bool fpe_trap::enable(void)
{
    // (unsupported platform.  leave fag set to false.
    // (using abortWithInsist to silence unused variable warning)
    if (abortWithInsist)
        fpeTrappingActive = false;
    else
        fpeTrappingActive = false;
    
    return fpeTrappingActive;
}

//---------------------------------------------------------------------------//
//! Disable trapping of floating point errors.
void fpe_trap::disable(void)
{
    return;
}

} // end namespace rtt_dsxx

#endif // FPETRAP_UNSUPPORTED

//---------------------------------------------------------------------------//
// end of fpe_trap.cc
//---------------------------------------------------------------------------//
