//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   fpe_trap/fpe_trap.hh
 * \author Rob Lowrie
 * \date   Thu Oct 13 16:36:09 2005
 * \brief  Contains functions in the fpe_trap namespace.
 * \note   Copyright (C) 2005-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef fpe_trap_hh
#define fpe_trap_hh

#include "ds++/config.h"

namespace rtt_fpe_trap
{

//---------------------------------------------------------------------------//
/*!
  \brief Enable trapping of floating-point exceptions.

  The floating-point exception behavior is platform dependent.  Nevertheless,
  the goal of this function is to turn on trapping for the following
  exceptions:

  - Division by zero.
  - Invalid operation; for example, the square-root of a negative
    number.
  - Overflow.

  If a floating-point exception is detected, ds++'s Insist is called; that is,
  a C++ exception is thrown.

  Typically, an application calls this function once, before any
  floating-point operations are performed.  Note that all program
  functionality then traps floating-point exceptions, including in libraries.
  Currently, there is no way to turn trapping off.

  \return \b true if trapping is enabled, \b false otherwise.
  A \b false return value is typically because the platform is not
  supported.
*/
DLL_PUBLIC bool enable_fpe();

} // end namespace rtt_fpe_trap

#endif // fpe_trap_hh

//---------------------------------------------------------------------------//
// end of fpe_trap/fpe_trap.hh
//---------------------------------------------------------------------------//
