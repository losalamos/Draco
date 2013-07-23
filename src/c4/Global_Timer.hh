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

#include <map>
#include <set>
#include <string>

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
 * All Global_Timers have a unique name assigned via the destructor. They can
 * be can be enabled or disabled as a whole by setting or unsetting a single
 * global variable, or they can be activated selectively by name. When
 * enabled, they automatically generate a report to cout when they go out of
 * scope. The class is carefully designed to avoid order of initialization or
 * destruction issues. Since instances of this class are normally declared at
 * global scope, the timing reports appear at the end of a calculation.
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

    bool active_;

    static bool global_active_; // All Global_Timers are active

    struct timer_entry
    {
        bool is_active;
        Global_Timer *timer;

        timer_entry() : is_active(false), timer(NULL) {}
    };

    static std::map<std::string, timer_entry> active_list_; // Selected Global_Timers are active
     
  public:
    
    Global_Timer(char const *name); //! default constructor

    char const *name() const { return name_; }

    bool is_active() const { return active_; }

    void set_activity(bool active) { active_ = active; }

    static bool is_global_active() { return global_active_; }

    static void set_global_activity(bool active);

    static void set_global_activity(std::set<std::string> const &timer_list);

    ~Global_Timer();
};

} // end namespace rtt_c4


#endif // __c4_Global_Timer_hh__

//---------------------------------------------------------------------------//
//                              end of c4/Global_Timer.hh
//---------------------------------------------------------------------------//
