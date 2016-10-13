//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/ParallelUnitTest.cc
 * \author Kelly Thompson
 * \date   Thu Jun  1 17:15:05 2006
 * \brief  Implementation file for encapsulation of Draco parallel unit tests.
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/XGetopt.hh"
#include <iostream>
#include <sstream>
#include <string>

namespace rtt_c4 {
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for ParallelUnitTest
 * \arg argc The number of command line arguments
 * \arg argv A list of strings containg the command line arguments
 * \arg release_ A function pointer to this package's release function.
 * \arg out_ A user specified iostream that defaults to std::cout.
 * \exception rtt_dsxx::assertion An exception with the message "Success" will
 * be thrown if \c --version is found in the argument list.
 *
 * The constructor initializes the prallel communicator (MPI) and then
 * initializes the base class UnitTest by setting numPasses and numFails to
 * zero.  It also prints a message that declares this to be a scalar unit test
 * and provides the unit test name.
 */
ParallelUnitTest::ParallelUnitTest(int &argc, char **&argv,
                                   string_fp_void release_, std::ostream &out_)
    : UnitTest(argc, argv, release_, out_) {
  using std::string;

  initialize(argc, argv);

  Require(argc > 0);
  Require(release != NULL);

  // header

  if (node() == 0)
    out << "\n============================================="
        << "\n=== Parallel Unit Test: " << testName
        << "\n=== Number of Processors: " << nodes()
        << "\n=============================================\n"
        << std::endl;

  // version tag

  if (node() == 0)
    out << testName << ": version " << release() << "\n" << std::endl;

  // Register and process command line arguments:
  rtt_dsxx::XGetopt::csmap long_options;
  long_options['v'] = "version";
  std::map<char, std::string> help_strings;
  help_strings['v'] = "print version information and exit.";
  rtt_dsxx::XGetopt program_options(argc, argv, long_options, help_strings);

  int c(0);
  while ((c = program_options()) != -1) {
    switch (c) {
    case 'v': // --version
      finalize();
      throw rtt_dsxx::assertion(string("Success"));
      break;

    default:
      return; // nothing to do.
    }
  }

  Ensure(numPasses == 0);
  Ensure(numFails == 0);
  Ensure(testName.length() > 0);

  return;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 * The destructor provides a final status report before it calls MPI_Finalize
 * and exits.
 *
 * Use global_sum to ensure that we print FAIL if any tests on any processor
 * fail.
 */
ParallelUnitTest::~ParallelUnitTest() {
  global_sum(numPasses);
  global_sum(numFails);
  if (node() == 0)
    out << resultMessage() << std::endl;
  // global_barrier();
  // out << resultMessage() << std::endl;
  global_barrier();
  finalize();
  return;
}

//---------------------------------------------------------------------------//
//! Print a summary of the pass/fail status of ParallelUnitTest.
void ParallelUnitTest::status() {
  { // Provide some space before the report -- but keep all the processors
    // in sync.  [KT: 2011/06/20 - Actually, ParallelUnitTest should only
    // have a barrier on the destructor.  Otherwise, we can find ourselves
    // in a race condition between this function and the destructor (in the
    // case of an exception).]
    if (node() == 0)
      out << std::endl;
  }
  {
    out << "Done testing " << testName << " on node " << node() << "."
        << std::endl;
  }
  return;
}

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
// end of ParallelUnitTest.cc
//---------------------------------------------------------------------------//
