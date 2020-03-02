//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   c4/Global_Timer.hh
 * \author Kent G. Budge
 * \date   Mon Mar 25 17:35:07 2002
 * \brief  Define class Global_Timer, a POSIX standard timer.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#ifndef rtt_c4_Global_Timer_hh
#define rtt_c4_Global_Timer_hh

#include "Timer.hh"
#include <map>
#include <set>

namespace rtt_c4 {

//============================================================================//
/*!
 * \class Global_Timer
 *
 * \brief POSIX standard timer.
 *
 * The Global_Timer class is based on the Timer class, but adds additional
 * features that make it more convenient for timing sections of code without
 * being tied to specific objects that use those sections of code.
 *
 * All Global_Timers have a unique name assigned via the constructor. They can
 * be enabled or disabled as a whole by setting or unsetting a single global
 * variable, or they can be activated selectively by name. The timings can be
 * reset manually, and reports for all active timers can be generated with a
 * single static function call.
 *
 * Global_Timers are only active on processor 0.
 */
//============================================================================//

class Global_Timer : public Timer {
private:
  char const *name_; // name assigned by client to this timer, to distinguish
                     // its output from that of any other timers.

  bool active_; // This timer is active. This does not mean it is currently
                // accumulating timing statistics, but only that it is flagged
                // to do so when start() is called. If not active, a call to
                // start() is ignored.

  //! All Global_Timers are active
  DLL_PUBLIC_c4 static bool global_active_;

  struct timer_entry {
    bool is_active; // permits activation of timers not yet constructed.
    Global_Timer *timer;

    timer_entry() : is_active(false), timer(nullptr) {}
  };

  typedef std::map<std::string, timer_entry> active_list_type;

  //! Selected Global_Timers are active
  static active_list_type active_list_;

  //! Disable default construction
  Global_Timer() = delete;

  //! Disable copy/move construction
  Global_Timer(Global_Timer const &rhs) = delete;
  Global_Timer(Global_Timer &&rhs) = delete;

  // Disable copy/move assignment
  Global_Timer operator=(Global_Timer const &rhs) = delete;
  Global_Timer operator=(Global_Timer &&rhs) = delete;

public:
  // Constructors & Destructors

  explicit Global_Timer(char const *name); //! default constructor
  virtual ~Global_Timer(){};               //! default destructor

  // Accessors

  char const *name() const { return name_; }

  bool is_active() const { return active_ || global_active_; }

  // Manipulators

  void set_activity(bool active) { active_ = active; }

  void start() {
    if (active_ || global_active_)
      Timer::start();
  }

  void stop() {
    if (active_ || global_active_)
      Timer::stop();
  }

  // Statics

  /*!
   *  \bug no documentation!
   */
  static void set_global_activity(bool active);

  static void set_selected_activity(std::set<std::string> const &timer_list,
                                    bool active);

  static void reset_all();

  static void report_all(std::ostream &);
};

} // end namespace rtt_c4

#endif // rtt_c4_Global_Timer_hh

//----------------------------------------------------------------------------//
// end of c4/Global_Timer.hh
//----------------------------------------------------------------------------//
