//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/windows_x86.cc
 * \author Rob Lowrie
 * \date   Thu Oct 13 16:52:05 2005
 * \brief  windows implementation of fpe_trap functions.
 * \note   Copyright (C) 2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 * \sa Microsoft Visual C++ Floating-Point Optimizations, 
 *     section "Floating-Point Exceptions as C++ Exceptions"
 */
//---------------------------------------------------------------------------//
// $Id: WINDOWS_x86.cc 5571 2010-12-10 20:07:38Z kellyt $
//---------------------------------------------------------------------------//

#include "fpe_trap/config.h"

#ifdef FPETRAP_WINDOWS_X86

#include <Windows.h>  // defines STATUS_FLOAT_...
#include <float.h>    // defines _controlfp_s
#include <iostream>
#include <string>
#include <ds++/Assert.hh>

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

namespace rtt_fpe_trap
{

// ----------------------------------------------------------------------------
// - http://stackoverflow.com/questions/2769814/how-do-i-use-try-catch-to-catch-floating-point-errors
// - See MSDN articles on fenv_access and _controlfp_s examples.
// ----------------------------------------------------------------------------
DLL_PUBLIC bool enable_fpe()
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
   
   _set_se_translator(trans_func);
   
   return true;
}

} // end namespace rtt_shared_lib

#endif // FPETRAP_WINDOWS_X86

//---------------------------------------------------------------------------//
// end of windows_x86.cc
//---------------------------------------------------------------------------//
