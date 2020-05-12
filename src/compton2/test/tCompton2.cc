//-----------------------------------*-C++-*----------------------------------//
/*!
 * \file   compton2/test/tCompton2.cc
 * \author Andrew Till
 * \date   11 May 2020
 * \brief  Implementation file for tCompton2
 * \note   Copyright (C) 2017-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "compton2/Compton2.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace rtt_compton2_test {

using rtt_dsxx::soft_equiv;

//----------------------------------------------------------------------------//
// TESTS
//----------------------------------------------------------------------------//

//!  Simple test
void test(rtt_dsxx::UnitTest &ut) {
  // Start the test.

  if (!soft_equiv(1.0, 1.0))
    FAILMSG("Very bad");

  if (ut.numFails == 0) {
    PASSMSG("Successly did nothing.");
  } else {
    FAILMSG("Did not successfully do nothing.");
  }
}

} // namespace rtt_compton2_test

//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    // >>> UNIT TESTS
    rtt_compton2_test::test(ut);
  }
  UT_EPILOG(ut);
}

//----------------------------------------------------------------------------//
// End of test/tCompton2.cc
//----------------------------------------------------------------------------//
