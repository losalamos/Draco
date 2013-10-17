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
  \class fpe_trap

  \brief Enable trapping of floating-point exceptions.
  \param[in] abortWithInsist toggle the abort mode default is true to use
                             ds++'s Insist macro.
                             
  The floating-point exception behavior is platform dependent.  Nevertheless,
  the goal of this class is to turn on trapping for the following exceptions:

  - Division by zero.
  - Invalid operation; for example, the square-root of a negative number.
  - Overflow.

  If a floating-point exception is detected, the code will abort using a mode
  triggered by the value of abortWithInsist.
  - If true, ds++'s Insist is called; that is, a C++ exception is thrown.
  - If false, the default mechanism defined by the compiler will be used.  For
    most modern compilers this results in a stack trace.

  Typically, an application calls this function once, before any
  floating-point operations are performed (e.g.: \c
  wedgehog/Function_Interfaces.cc).  Note that all program functionality then
  traps floating-point exceptions, including in libraries.  Currently, there
  is no way to turn trapping off once it has been turned on.

  Note: By Draco coding convention, fpe_traps are enabled when
  DRACO_DIAGNOSTIC && 4 == true.
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


    //-----------------------------------------------------------------------------------//
    /*! 
     * \brief Enable trapping fpe signals.
     * \return \b true if trapping is enabled, \b false otherwise.
     *         A \b false return value is typically because the platform is not 
     *         supported.
     */
    bool enable( void );
    //! Disable trapping fpe signals.
    void disable( void );
    //! Query if trapping of fpe signals is active.
    bool active(void) const { return fpeTrappingActive; }

  private:
    bool fpeTrappingActive;
    bool const abortWithInsist;
};

} // end namespace rtt_fpe_trap

#endif // fpe_trap_hh

//---------------------------------------------------------------------------//
// end of fpe_trap/fpe_trap.hh
//---------------------------------------------------------------------------//
