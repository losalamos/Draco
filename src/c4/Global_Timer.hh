//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Global_Timer.hh
 * \author Kent G. Budge
 * \date   Mon Mar 25 17:35:07 2002
 * \brief  Define class Global_Timer, a POSIX standard timer.
 * \note   Copyright (C) 2002-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: Timer.hh 7075 2013-04-01 22:48:15Z kellyt $
//---------------------------------------------------------------------------//

#ifndef __c4_Global_Timer_hh__
#define __c4_Global_Timer_hh__

#include "Timer.hh"

namespace rtt_c4
{

//===========================================================================//
/*!
 * \class Global_Timer
 *
 * \brief POSIX standard timer.
 *
 * The Global_Timer class is based on the Timer class, but adds additional
 * features that make it more convenient for timing sections of code without
 * being tied to specific objects that use those sections of code.
 *
 * The chief of these features is that all Global_Timers are enabled or
 * disabled by a single global variable, and, when enabled, they
 * automatically generate a report to cout when they go out of scope. Since
 * instances of this class are normally declared at global scope, this report
 * appears at the end of a calculation.
 *
 * Global_Timers are only active on processor 0.
 */
//===========================================================================//

class DLL_PUBLIC Global_Timer : public Timer
{
  private:

    char const *name_; // name assigned by client to this timer, to
                       // distinguish its output from that of any other
                       // timers.

    static bool active_; // Global_Timers are active
     
  public:
    
    Global_Timer(char const *name); //! default constructor

    char const *name() const { return name_; }

    static bool is_active() { return active_; }

    static void set_activity(bool active);

    ~Global_Timer();
};

} // end namespace rtt_c4


#endif // __c4_Global_Timer_hh__

//---------------------------------------------------------------------------//
//                              end of c4/Global_Timer.hh
//---------------------------------------------------------------------------//
