//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/fpe_trap.hh
 * \author Rob Lowrie
 * \date   Thu Oct 13 16:36:09 2005
 * \brief  Contains functions in the fpe_trap namespace.
 * \note   Copyright (C) 2005-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef fpe_trap_hh
#define fpe_trap_hh

#include "ds++/config.h"
#include <string>
#include <iostream>

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
/*!
 * \class fpe_trap
 *
 * \brief Enable trapping of floating-point exceptions.
 * \param[in] abortWithInsist toggle the abort mode default is true to use
 *                            ds++'s Insist macro.
 *                            
 * The floating-point exception behavior is platform dependent.  Nevertheless,
 * the goal of this class is to turn on trapping for the following exceptions:
 *
 * - Division by zero.
 * - Invalid operation; for example, the square-root of a negative number.
 * - Overflow.
 *
 * If a floating-point exception is detected, the code will abort using a mode
 * triggered by the value of abortWithInsist.
 * - If true, ds++'s Insist is called; that is, a C++ exception is thrown.
 * - If false, the default mechanism defined by the compiler will be used.  For
 *   most modern compilers this results in a stack trace.
 *
 * Typically, an application calls this function once, before any
 * floating-point operations are performed (e.g.: \c
 * wedgehog/Function_Interfaces.cc).  Note that all program functionality then
 * traps floating-point exceptions, including in libraries.  Currently, there
 * is no way to turn trapping off once it has been turned on.
 *
 * Note: By Draco coding convention, fpe_traps are enabled when
 * DRACO_DIAGNOSTIC && 4 == true.
 *
 * Useful links:
 * - http://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes
 */
class DLL_PUBLIC fpe_trap
{
  public:
    //! constructor
    fpe_trap( bool const abortWithInsist_in = true)
        : fpeTrappingActive(false),
          abortWithInsist( abortWithInsist_in )
    {/* emtpy */};
    ~fpe_trap(void) {/* empty */};

    //! Enable trapping of fpe signals.
    bool enable( void );
    //! Disable trapping of fpe signals.
    void disable( void );
    //! Query if trapping of fpe signals is active.
    bool active(void) const { return fpeTrappingActive; }

  private:
    bool fpeTrappingActive;
    bool abortWithInsist;
};

} // end namespace rtt_dsxx

#endif // fpe_trap_hh

//---------------------------------------------------------------------------//
// end of ds++/fpe_trap.hh
//---------------------------------------------------------------------------//
