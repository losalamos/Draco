//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/tstRnd_Control_Inline.cc
 * \author Paul Henning
 * \brief  Rnd_Control test.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "rng/Rnd_Control_Inline.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_rng;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_control(UnitTest &ut) {
  // Create two controllers with different seeds.
  uint32_t seed1 = 2452423;
  uint32_t seed2 = 9182736;
  Rnd_Control control1(seed1);
  Rnd_Control control2(seed2);

  if (control1.get_max_streams() != std::numeric_limits<uint64_t>::max())
    ITFAILS;
  if (control2.get_max_streams() != std::numeric_limits<uint64_t>::max())
    ITFAILS;
  if (control1.get_seed() != seed1)
    ITFAILS;
  if (control2.get_seed() != seed2)
    ITFAILS;
  if (control1.get_num() != 0)
    ITFAILS;
  if (control2.get_num() != 0)
    ITFAILS;

  // Create a third controller with the same seed as controller1, but starting
  // with a different stream number.
  uint64_t streamnum = 2000;
  Rnd_Control control3(seed1, streamnum);

  // Initialize some generators.
  Counter_RNG rng1, rng2, rng3;
  uint64_t numiter = 1000;
  for (uint64_t i = 0; i < numiter; ++i) {
    control1.initialize(rng1);
    control2.initialize(rng2);
    control3.initialize(rng3);

    // Both rng1 and rng2 should be on stream number i.  control1 and control2
    // should be on stream number i+1.  rng3 should be on stream number
    // i+streamnum.  control3 should be on stream number i+streamnum+1.  Other
    // controller state should not have changed.  None of the generators should
    // match each other.
    if (rng1.get_num() != i)
      ITFAILS;
    if (rng2.get_num() != i)
      ITFAILS;
    if (rng3.get_num() != i + streamnum)
      ITFAILS;

    if (rng1 == rng2)
      ITFAILS;
    if (rng1 == rng3)
      ITFAILS;
    if (rng2 == rng3)
      ITFAILS;

    if (control1.get_num() != i + 1)
      ITFAILS;
    if (control2.get_num() != i + 1)
      ITFAILS;
    if (control3.get_num() != i + streamnum + 1)
      ITFAILS;

    if (control1.get_max_streams() != std::numeric_limits<uint64_t>::max())
      ITFAILS;
    if (control2.get_max_streams() != std::numeric_limits<uint64_t>::max())
      ITFAILS;
    if (control3.get_max_streams() != std::numeric_limits<uint64_t>::max())
      ITFAILS;
    if (control1.get_seed() != seed1)
      ITFAILS;
    if (control2.get_seed() != seed2)
      ITFAILS;
    if (control3.get_seed() != seed1)
      ITFAILS;
  }

  // Create another controller with the same seed and control number as the
  // original controller.
  Rnd_Control control4(seed1);

  if (control4.get_num() != 0)
    ITFAILS;

  // Set control4's next stream number manually.
  control4.set_num(numiter - 1);

  if (control4.get_num() != numiter - 1)
    ITFAILS;

  // Initialize a generator.
  Counter_RNG rng4;
  control4.initialize(rng4);

  // rng4 should match rng1 but not rng2 or rng3.
  if (rng4.get_num() != numiter - 1)
    ITFAILS;
  if (rng4 != rng1)
    ITFAILS;
  if (rng4 == rng2)
    ITFAILS;
  if (rng4 == rng3)
    ITFAILS;
  if (control4.get_num() != numiter)
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("test_control passed");
}

//---------------------------------------------------------------------------//
void test_exceptions(UnitTest &ut) {
// 1. Only test exceptions if DbC is enabled.
// 2. However, do not run these tests if no-throw DbC is enabled (DBC & 8)
#ifdef REQUIRE_ON
#if !(DBC & 8)
  // Try to create a controller that allows 0 streams.
  bool caught = false;
  try {
    Rnd_Control control(0, 0, 0);
  } catch (rtt_dsxx::assertion &err) {
    cout << "Good, caught assertion: " << err.what() << endl;
    caught = true;
  }
  if (!caught)
    ITFAILS;

  // Try to create a controller with an initial stream number greater than its
  // maximum number of streams.
  caught = false;
  try {
    Rnd_Control control(0, 1001, 1000);
  } catch (rtt_dsxx::assertion &err) {
    cout << "Good, caught assertion: " << err.what() << endl;
    caught = true;
  }
  if (!caught)
    ITFAILS;

  // Create a controller.
  Rnd_Control control(0, 0);

  // Try to set the stream number to std::numeric_limits<uint64_t>::max().
  caught = false;
  try {
    control.set_num(std::numeric_limits<uint64_t>::max());
  } catch (rtt_dsxx::assertion &err) {
    cout << "Good, caught assertion: " << err.what() << endl;
    caught = true;
  }
  if (!caught)
    ITFAILS;

  // Set the stream number to std::numeric_limits<uint64_t>::max() - 1, then try
  // to initialize two generators.  One should succeed.
  control.set_num(std::numeric_limits<uint64_t>::max() - 1);
  caught = false;
  uint64_t num_rngs = 0;
  try {
    for (num_rngs = 0; num_rngs < 2; ++num_rngs) {
      Counter_RNG rng;
      control.initialize(rng);

      if (rng.get_num() != std::numeric_limits<uint64_t>::max() - 1)
        ITFAILS;
    }
  } catch (rtt_dsxx::assertion &err) {
    cout << "Good, caught assertion: " << err.what() << endl;
    caught = true;
  }
  if (!caught)
    ITFAILS;
  if (num_rngs != 1)
    ITFAILS;

  // Create a controller that allows 10 streams.
  Rnd_Control control2(0, 0, 10);

  if (control2.get_max_streams() != 10)
    ITFAILS;

  // Try to create a generator on stream 10.
  caught = false;
  try {
    Counter_RNG rng;
    control2.initialize(10, rng);
  } catch (rtt_dsxx::assertion &err) {
    cout << "Good, caught assertion: " << err.what() << endl;
    caught = true;
  }
  if (!caught)
    ITFAILS;

  // Try to create 11 generators.  The first 10 should succeed.
  num_rngs = 0;
  caught = false;
  try {
    for (num_rngs = 0; num_rngs < 11; ++num_rngs) {
      Counter_RNG rng;
      control2.initialize(rng);

      if (rng.get_num() != num_rngs)
        ITFAILS;
    }
  } catch (rtt_dsxx::assertion &err) {
    cout << "Good, caught assertion: " << err.what() << endl;
    caught = true;
  }
  if (!caught)
    ITFAILS;
  if (num_rngs != 10)
    ITFAILS;
#endif
#endif

  if (ut.numFails == 0)
    PASSMSG("test_exceptions passed");
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    test_control(ut);
    test_exceptions(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstRnd_Control.cc
//---------------------------------------------------------------------------//
