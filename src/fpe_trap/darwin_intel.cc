//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/darwin_intel.cc
 * \author Rob Lowrie
 * \date   Sun Jul 19 21:43:45 2009
 * \brief  OS/X Intel implementation of fpe_trap functions.
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

#include <fpe_trap/config.h>

#ifdef FPETRAP_DARWIN_INTEL

#include <iostream>
#include <string>
#include <ds++/Assert.hh>

#include <xmmintrin.h>

namespace rtt_fpe_trap
{

//---------------------------------------------------------------------------//
//! Enable trapping of floating point errors.
bool fpe_trap::enable(void)
{
    _mm_setcsr( _MM_MASK_MASK &~
                (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );    

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
// end of darwin_intel.cc
//---------------------------------------------------------------------------//
