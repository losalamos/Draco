//----------------------------------*-C++-*----------------------------------//
/*! 
 * \file   timestep/Release.hh
 * \author John M. McGhee
 * \date   Fri Aug 27 10:33:26 1999
 * \brief  Header file for timestep library release function.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __timestep_Release_hh__
#define __timestep_Release_hh__

//===========================================================================//

#include <string>

/*!
 * \brief RTT time step namespace.
 *
 * Provides namespace protection for the Draco (RTT) time step 
 * control utilities.
 *\sa The ts_manager and ts_advisor classes provide most of the
 *    functionality of the namespace. The \ref timestep_overview page presents
 *    a summary of the capabilities provided within the namespace.
 */
namespace rtt_timestep 
{
/*!
 * \brief  Gets the release number for the timestep package. 
 * \return release number as a string in the form "timestep-\#_\#_\#"
 */
    const std::string release();
}

#endif                          // __timestep_Release_hh__


/*!
 * \page timestep_overview Overview of the Draco Time Step Manager
 *
 * \version 1_0_1
 *
 * <h3> Introduction </h3>
 * The classes contained in the rtt_timestep name space are designed
 * to provide automatic control of the time step size in a time dependent
 * numerical simulation. In a typical simulation, several physics
 * packages may be repeatedly called from a host code to advance the
 * overall system state. 
 * A time step manager class (rtt_timestep::ts_manager) is provided 
 * to contain any number of
 * time step advisors (rtt_timestep::ts_advisor) from any package 
 * in the simulation. Each advisor
 * provides a recommended time step based on some physical parameter
 * specific to the package. After collecting information from each
 * advisor, the time step manager considers all the inputs and provides
 * an overall recommendation for the size of the time step for the
 * next time cycle.  Capabilities are provided to:
 * <ul>
 *  <li> limit the min and max timestep,
 *  <li> force the time step to a user input value,
 *  <li> display informational advisors,
 *  <li> print a variety of summaries,
 *  <li> limit the rate of change in the time step
 *  <li> control the time step based on estimated rate of change in
 *       a user selected physical parameter,
 *  <li> attempt to hit a "target" time. 
 * </ul> 
 *
 * <h3> Intended Usage </h3>
 * The intended usage is as follows:
 *
 * 1) Host-code instantiates the desired global (i.e. package
 *    independent) time-step advisor(s).
 *
 * 2) Host-code instantiates a ts_manager, and loads references to
 *    advisor(s) created in 1) above into the manager list.
 *
 * 3) Each package in the host code, at package construction, 
 *    is passed in the manager, and instantiates any 
 *    desired time-step advisors and loads references 
 *    to the new advisors into the time-step manager.
 *
 * ***Begin loop over time ***
 *
 *    4) Host-code responds to any user query about advisor(s), and
 *       sets advisor control parameters, based on user interactive
 *       input or some pre-determined script.
 *
 *    5) Host-code sets cycle data of manager (dt, cycle#, time).
 *
 *    6) Each package  updates, if required, the recommended time
 *       step for each advisor "owned" by that package.
 *
 *    7) Host-code updates recommended time step(s) of global 
 *       advisor(s) not owned by any particular package, if required.
 *
 *    8) Manager evaluates recommendations from each advisor on the
 *       list, calculates new time-step, and prints results
 *       if desired.
 *
 * *** End loop over time ***
 *
 * 9)  At package destruction, each package removes the advisors
 *     created in 3) above from the manager.
 *
 * 10) Host-code removes global advisors from the manager list.
 *
 * 11) Host-code destructs manager.
 *
 * <h3> Other Draco Packages </h3>
 * The time step manager utility uses the Draco ds++ services library, 
 * and the Draco C4 communication library.
 */

/*!
 * \example timestep/test/main.cc
 * The following code provides an example of how to use the timestep manager
 * utility. It includes a dummy_package for use with the manager. Also
 * included isis a sample output from the test_timestep example.
 * It contains representative output from most of the printing and summary
 * I/O utilities.
 * \include timestep/test/test_timestep.hh
 * \include timestep/test/test_timestep.cc
 * \include timestep/test/dummy_package.hh
 * \include timestep/test/test_timestep_pt.cc
 * \include timestep/test/dummy_package.cc
 * \include timestep/test/test_timestep.out
 */

//---------------------------------------------------------------------------//
//                              end of timestep/Release.hh
//---------------------------------------------------------------------------//
