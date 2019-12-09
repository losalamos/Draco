//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/Timing.hh
 * \author T.M. Kelly, Thomas M. Evans
 * \date   Tue Dec 13 10:44:29 2005
 * \brief  Timing class and macros definition.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef diagnostics_Timing_hh
#define diagnostics_Timing_hh

#include "diagnostics/config.h"
#include "ds++/config.h"
#include <map>
#include <string>
#include <vector>

namespace rtt_diagnostics {

//===========================================================================//
/*!
 * \class Timing_Diagnostics
 * \brief Class to hold timing results for diagnostic output.
 *
 * This class provides a simple interface to store timing data during a
 * simulation.  Generally, one adds a timer label and value using the
 * update_timer() function:
 * \code
 *   // start a timer
 *   ...
 *   // do work
 *   ...
 *   // stop timer
 *   rtt_diagnostics::Timing_Diagnostics::update_timer("Solver", time);
 * \endcode
 * There is no need to "add" the timer entry for "Solver" before an update.
 * If the key "Solver" does not exist it is added with value 0.0 before
 * applying the value.
 *
 * The easiest way to use this class is through the TIMER macros.
 */
/*!
 * \example diagnostics/test/tstTiming.cc
 */
//===========================================================================//

class DLL_PUBLIC_diagnostics Timing_Diagnostics {
public:
  // Useful typedef.
  typedef std::vector<std::string> Vec_Keys;

private:
  // >>> PRIVATE DATA MEMBERS

  //! Map of timers.
  static std::map<std::string, double> timers;

public:
  // >>> FUNCTIONAL INTERFACE

  // Add a value to the timer with name key.
  static void update_timer(const std::string &key, double value);

  //! Get a timer's value.  Adds timer with name key to map.
  static double timer_value(const std::string &k) { return timers[k]; }

  //! Get number of timers in map.
  static size_t num_timers() { return timers.size(); }

  // Return a vector of timer keys.
  static Vec_Keys timer_keys();

  // Reset a timer.
  static void reset_timer(const std::string &key);

  // Reset all timers from the map of timers.
  static void reset_timers();

  // Delete a timer from the map of timers.
  static void delete_timer(const std::string &key);

  // Delete all timers from the map of timers.
  static void delete_timers();

private:
  // >>> IMPLEMENTATION

  // This class is never constructed.
  Timing_Diagnostics();

  // This class is also never destructed.
  ~Timing_Diagnostics();
};

} // end namespace rtt_diagnostics

//---------------------------------------------------------------------------//
/*!
 * \page diagnostics_timing Macros for timing
 *
 * Four macros are defined here; these macros insert timer calls into code if
 * a global definition, DRACO_TIMING, is greater then 0.  The default value is
 * to set DRACO_TIMING == 0. They use the draco rtt_c4::Timer and
 * rtt_diagnostics::Timing_Diagnostics classes.
 *
 * The build system sets DRACO_TIMING through the configure option \c
 * --with-clubimc-timing.  The following settings apply:
 * - 0 turns off all TIMER macros
 * - 1 turns on TIMER, TIMER_START, TIMER_STOP, and TIMER_RECORD
 * - 2 turns on all TIMER macros (include TIMER_REPORT)
 * .
 * The default is 0.
 *
 * In code,
 * \code
 * #include "diagnostics/Timing.hh"
 *
 * TIMER( foo);
 * TIMER_START( foo);
 * // ...
 * // code interval to time
 * // ...
 * TIMER_STOP( foo);
 * TIMER_RECORD( "Snippet", foo);
 * TIMER_REPORT( foo, std::cout, "interval 42");
 * \endcode
 * would produce the following output:
 * \verbatim
 *   somefile.cc ###: interval 42 elapsed wall_clock: ## seconds; elapsed user_time: ## seconds; elapsed sys_time ## seconds.\n
 * \endverbatim
 * The key "Snippet" can be used to access the stored time through the
 * Timer_Diagnostics class:
 * \code
 * #ifdef DRACO_TIMING_ON
 *   vector<string> keys = Timing_Diagnostics::timer_keys();
 *   for (int i = 0; i < keys.size(); i++)
 *   {
 *       cout << keys[i] << "\t" << Timing_Diagnostics::timer_value(keys[i])
 *            << endl;
 *   }
     Timing_Diagnostics::reset_timers();
 * #endif
 * \endcode
 */

/*!
 * \def TIMER( timer_name)
 *
 * If DRACO_TIMING_ON is defined, TIMER( timer_name) expands to:
 * \code
 *     rtt_c4::Timer timer_name
 * \endcode
 * Otherwise it is empty.
 */

/*!
 * \def TIMER_START( timer_name)
 *
 * If DRACO_TIMING > 0 TIMER_START( timer_name) expands to:
 * \code
 *     timer_name.start()
 * \endcode
 * Otherwise it is empty.
 */

/*!
 * \def TIMER_STOP( timer_name)
 *
 * If DRACO_TIMING_ON > 0, TIMER_STOP( timer_name) expands to:
 * \code
 *     timer_name.stop()
 * \endcode
 * Otherwise it is empty.
 */

/*!
 * \def TIMER_RECORD( name, timer)
 *
 * If DRACO_TIMING_ON > 0, TIMER_RECORD( name, timer) expands to:
 * \code
 *     rtt_diagnostics::Timing_Diagnostics::update_timer(name, timer.wall_clock())
 * \endcode
 * Otherwise it is empty.
 */

/*!
 * \def TIMER_REPORT( timer_name, ostream, comment)
 *
 * If DRACO_TIMING > 1, TIMER_REPORT( timer_name, ostream,
 * comment) expands to:
 * \code
 *     ostream << __FILE__ << " " << __LINE__ << ": " << comment      \
 *             << " elapsed wall_clock: " << timer.wall_clock()       \
 *             << " seconds; elapsed user_time: " << timer.user_cpu() \
 *             << " seconds; elapsed sys_time: " << timer.system_cpu()\
 *             << " seconds.\n" << flush
 * \endcode
 * Otherwise it is empty. The flush ensures that regression tests
 * continue to pass (otherwise, in parallel runs, output may arrive "out of
 * order" and trample the output that the regression tests look for).
 */
//---------------------------------------------------------------------------//

#if !defined(DRACO_TIMING)
#define DRACO_TIMING 0
#endif

//---------------------------------------------------------------------------//
/*
 * All timing operations are inactive.
 */
#if DRACO_TIMING == 0

#define TIMER(timer)

#define TIMER_START(timer)

#define TIMER_STOP(timer)

#define TIMER_RECORD(name, timer)

#define TIMER_REPORT(timer, ostream, comment)

#endif

//---------------------------------------------------------------------------//
/*
 * Turn on basic timing operations.
 */
#if DRACO_TIMING > 0

#include "c4/Timer.hh"

#define DRACO_TIMING_ON

#define TIMER(timer) rtt_c4::Timer timer

#define TIMER_START(timer) timer.start()

#define TIMER_STOP(timer) timer.stop()

#define TIMER_RECORD(name, timer)                                              \
  rtt_diagnostics::Timing_Diagnostics::update_timer(name, timer.wall_clock())

#endif

//---------------------------------------------------------------------------//
/*
 * Turn on timing report output.  This is an add-on option to the basic timing
 * operations.
 */
#if DRACO_TIMING > 1

#define TIMER_REPORT(timer, ostream, comment)                                  \
  ostream << __FILE__ << " " << __LINE__ << ": " << comment                    \
          << " elapsed wall_clock: " << timer.wall_clock()                     \
          << " seconds; elapsed user_time: " << timer.user_cpu()               \
          << " seconds; elapsed sys_time: " << timer.system_cpu()              \
          << " seconds.\n"                                                     \
          << std::flush

#else

#define TIMER_REPORT(timer, ostream, comment)

#endif

#endif // diagnostics_Timing_hh

//---------------------------------------------------------------------------//
//  end of diagnostics/Timing.hh
//---------------------------------------------------------------------------//
