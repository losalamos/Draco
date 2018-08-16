//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/tstCounter_RNG2.cc
 * \author Peter Ahrens
 * \date   Fri Aug 3 16:53:23 2012
 * \brief  Counter_RNG2 tests.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "rng/Counter_RNG2.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include <set>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_rng;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

/* I've included some test of getters and setters since they're not trivial. */
void test_set_get_seed(UnitTest &ut){
  uint32_t seed = 1;
  uint64_t const streamnum = 2;
  uint64_t const time_step = 349;
  Counter_RNG2 rng(seed, streamnum, time_step);

  rng.ran();
  rng.ran();
  rng.ran();

  rng.set_seed(42);
  FAIL_IF_NOT(42 == rng.get_seed());
  FAIL_IF_NOT(3 == rng.get_step_counter());
  FAIL_IF_NOT(time_step == rng.get_time_step());
  FAIL_IF_NOT(streamnum == rng.get_stream_number());
  FAIL_IF_NOT(0 == rng.get_spawn_id());

  if (ut.numFails == 0)
    PASSMSG("test_set_get_seed passed");
} // test_set_get_seed

void test_set_get_step_counter(UnitTest &ut){
  uint32_t const seed = 1;
  uint64_t const streamnum = 2;
  uint64_t const time_step = 349;
  Counter_RNG2 rng(seed, streamnum, time_step);

  rng.ran();
  rng.ran();
  rng.ran();

  rng.set_step_counter(42);
  FAIL_IF_NOT(seed == rng.get_seed());
  FAIL_IF_NOT(42 == rng.get_step_counter());
  FAIL_IF_NOT(time_step == rng.get_time_step());
  FAIL_IF_NOT(streamnum == rng.get_stream_number());
  FAIL_IF_NOT(0 == rng.get_spawn_id());

  if (ut.numFails == 0)
    PASSMSG("test_set_get_step_counter passed");
} // test_set_get_step_counter

void test_set_get_time_step(UnitTest &ut){
  uint32_t const seed = 1;
  uint64_t const streamnum = 2;
  uint64_t time_step = 349;
  Counter_RNG2 rng(seed, streamnum, time_step);

  rng.ran();
  rng.ran();
  rng.ran();

  rng.set_time_step(1000000);
  FAIL_IF_NOT(seed == rng.get_seed());
  FAIL_IF_NOT(3 == rng.get_step_counter());
  FAIL_IF_NOT(1000000 == rng.get_time_step());
  FAIL_IF_NOT(streamnum == rng.get_stream_number());
  FAIL_IF_NOT(0 == rng.get_spawn_id());

  if (ut.numFails == 0)
    PASSMSG("test_set_get_time_step passed");
} // test_set_get_time_step

void test_set_get_stream_number(UnitTest &ut){
  uint32_t const seed = 1;
  uint64_t streamnum = 2;
  uint64_t const time_step = 349;
  Counter_RNG2 rng(seed, streamnum, time_step);

  rng.ran();
  rng.ran();
  rng.ran();

  rng.set_stream_number(1000000);
  FAIL_IF_NOT(seed == rng.get_seed());
  FAIL_IF_NOT(3 == rng.get_step_counter());
  FAIL_IF_NOT(time_step == rng.get_time_step());
  FAIL_IF_NOT(1000000 == rng.get_stream_number());
  FAIL_IF_NOT(0 == rng.get_spawn_id());

  if (ut.numFails == 0)
    PASSMSG("test_set_get_stream_number passed");
} // test_set_get_stream_number

void test_set_get_spawn_id(UnitTest &ut){
  uint32_t const seed = 1;
  uint64_t const streamnum = 2;
  uint64_t const time_step = 349;
  Counter_RNG2 rng(seed, streamnum, time_step);

  rng.ran();
  rng.ran();
  rng.ran();

  rng.set_spawn_id(1000000);
  FAIL_IF_NOT(seed == rng.get_seed());
  FAIL_IF_NOT(3 == rng.get_step_counter());
  FAIL_IF_NOT(time_step == rng.get_time_step());
  FAIL_IF_NOT(streamnum == rng.get_stream_number());
  FAIL_IF_NOT(1000000 == rng.get_spawn_id());

  if (ut.numFails == 0)
    PASSMSG("test_set_get_spawn_id passed");
} // test_set_get_spawn_id

void test_increment_CBRNG2_step_counter(UnitTest &ut){
  uint64_t const ctr1 = 0xFFFFFFFF00000000;
  uint64_t const inc1 = increment_CBRNG2_step_counter(ctr1);
  FAIL_IF_NOT(inc1 == ctr1 + 1);
  uint64_t const ctr2 = 0xFFFFFFFFFFFFFFFF;
  uint64_t const inc2 = increment_CBRNG2_step_counter(ctr2);
  FAIL_IF_NOT(inc2 == ctr1);
  if (ut.numFails == 0)
    PASSMSG("test_increment_CBRNG2_step_counter passed");
}

void test_equality(UnitTest &ut) {
  // Create a Counter_RNG2 by specifying a seed and stream number.
  uint32_t seed = 1;
  uint64_t streamnum = 2;
  uint64_t time_step = 349;
  Counter_RNG2 rng(seed, streamnum, time_step);

  if (rng.get_stream_number() != streamnum)
    ITFAILS;
  if (rng.size() != CBRNG_DATA_SIZE)
    ITFAILS;
  if (rng.size_bytes() != CBRNG_DATA_SIZE * sizeof(uint64_t))
    ITFAILS;
  if (rng != rng)
    ITFAILS;

  // Create another Counter_RNG2 with a different seed.
  seed = 2;
  Counter_RNG2 rng2(seed, streamnum, time_step);

  // rng2's stream number should match rng's, but the two generators should not
  // be identical.
  if (rng2.get_stream_number() != streamnum)
    ITFAILS;
  if (rng2.get_stream_number() != rng.get_stream_number())
    ITFAILS;
  if (rng2.get_unique_num() == rng.get_unique_num())
    ITFAILS;
  if (rng2 == rng)
    ITFAILS;
  if (rng2 != rng2)
    ITFAILS;

  // Create another Counter_RNG2 with a different stream number.
  seed = 1;
  streamnum = 3;
  Counter_RNG2 rng3(seed, streamnum, time_step);

  // rng3 should be different from the previous two generators.
  if (rng3.get_stream_number() != streamnum)
    ITFAILS;
  if (rng3.get_unique_num() == rng.get_unique_num())
    ITFAILS;
  if (rng3.get_unique_num() == rng2.get_unique_num())
    ITFAILS;
  if (rng3 == rng)
    ITFAILS;
  if (rng3 == rng2)
    ITFAILS;
  if (rng3 != rng3)
    ITFAILS;

  // Create another Counter_RNG2 with the original seed and stream number.
  streamnum = 2;
  Counter_RNG2 rng4(seed, streamnum, time_step);

  // rng4 should be equal to rng but different from rng2 and rng3.
  if (rng4.get_stream_number() != streamnum)
    ITFAILS;
  if (rng4.get_unique_num() != rng.get_unique_num())
    ITFAILS;
  if (rng4.get_unique_num() == rng2.get_unique_num())
    ITFAILS;
  if (rng4.get_unique_num() == rng3.get_unique_num())
    ITFAILS;
  if (rng4 != rng)
    ITFAILS;
  if (rng4 == rng2)
    ITFAILS;
  if (rng4 == rng3)
    ITFAILS;
  if (rng4 != rng4)
    ITFAILS;

  // Create a Counter_RNG2 from a data array.
  vector<uint64_t> data(CBRNG_DATA_SIZE);
  data[0] = 1234;
  data[1] = 5678;
  data[2] = 9012;
  data[3] = 3456;
  Counter_RNG2 rng5(&data[0], &data[0] + CBRNG_DATA_SIZE);

  streamnum = data[1];
  if (rng5.get_stream_number() != streamnum)
    ITFAILS;
  if (rng5.get_unique_num() == rng.get_unique_num())
    ITFAILS;
  if (rng5.get_unique_num() == rng2.get_unique_num())
    ITFAILS;
  if (rng5.get_unique_num() == rng3.get_unique_num())
    ITFAILS;
  if (rng5.get_unique_num() == rng4.get_unique_num())
    ITFAILS;
  if (rng5 == rng)
    ITFAILS;
  if (rng5 == rng2)
    ITFAILS;
  if (rng5 == rng3)
    ITFAILS;
  if (rng5 == rng4)
    ITFAILS;
  if (rng5 != rng5)
    ITFAILS;

  // TK: I cut a test that was based on manually manipulating the RNG
  // state. Pls don't manually manipulate the RNG state as an array of ints.

// Try to create a Counter_RNG2 from a data array that's too short.
// 1. Only test exceptions if DbC is enabled.
// 2. However, do not run these tests if no-throw DbC is enabled (DBC & 8)
#ifdef REQUIRE_ON
#if !(DBC & 8)
  bool caught = false;
  try {
    Counter_RNG2 rng7(&data[0], &data[0] + CBRNG_DATA_SIZE - 1);
  } catch (rtt_dsxx::assertion &err) {
    cout << "Good, caught assertion: " << err.what() << endl;
    caught = true;
  }
  if (!caught)
    ITFAILS;
#endif
#endif

  if (ut.numFails == 0)
    PASSMSG("test_equality passed");
}

//---------------------------------------------------------------------------//
void test_stream(UnitTest &ut) {
  // Create two identical Counter_RNG2s.
  uint32_t seed = 0x12121212;
  uint64_t streamnum = 1234;
  uint64_t time_step = 456789;
  Counter_RNG2 rng(seed, streamnum, time_step);
  Counter_RNG2 rng2(seed, streamnum, time_step);

  if (rng != rng2)
    ITFAILS;

  // Generate a random double (and advance the stream) from rng.
  double x = rng.ran();

  // rng and rng2 should no longer match, but their stream numbers and unique
  // identifiers should be the same.
  if (rng == rng2)
    ITFAILS;
  if (rng.get_stream_number() != streamnum)
    ITFAILS;
  if (rng.get_stream_number() != rng2.get_stream_number())
    ITFAILS;
  if (rng.get_unique_num() != rng2.get_unique_num())
    ITFAILS;

  // Generate a random double (and advance the stream) from rng2.
  double y = rng2.ran();

  // Now rng and rng2 should match again, and the two generated doubles should
  // be identical.
  if (rng != rng2)
    ITFAILS;
  if (!soft_equiv(x, y))
    ITFAILS;

  // Generate another random double from rng.
  double z = rng.ran();

  // Now they should differ again.
  if (rng == rng2)
    ITFAILS;
  if (rng.get_stream_number() != streamnum)
    ITFAILS;
  if (rng.get_stream_number() != rng2.get_stream_number())
    ITFAILS;
  if (rng.get_unique_num() != rng2.get_unique_num())
    ITFAILS;
  if (soft_equiv(x, z))
    ITFAILS;

  // Create a Counter_RNG2 from a data array.
  vector<uint64_t> data(CBRNG_DATA_SIZE);
  data[0] = 0;
  data[1] = streamnum;
  data[2] = seed;
  data[3] = time_step;
  Counter_RNG2 rng3(&data[0], &data[0] + CBRNG_DATA_SIZE);

  // Initially, rng3 should exactly match neither rng nor rng2, but all three
  // should have the same stream number and "unique" identifier.
  if (!std::equal(rng3.begin(), rng3.end(), data.begin()))
    ITFAILS;
  if (rng3 == rng)
    ITFAILS;
  if (rng3 == rng2)
    ITFAILS;
  if (rng3.get_stream_number() != streamnum)
    ITFAILS;
  if (rng3.get_unique_num() != rng.get_unique_num())
    ITFAILS;
  if (rng3.get_unique_num() != rng2.get_unique_num())
    ITFAILS;

  // Generate a random double from rng3; it should match rng2 but not data
  // afterward.
  double w = rng3.ran();
  if (rng3 == rng)
    ITFAILS;
  if (rng3 != rng2)
    ITFAILS;
  if (std::equal(rng3.begin(), rng3.end(), data.begin()))
    ITFAILS;
  if (rng3.get_stream_number() != streamnum)
    ITFAILS;
  if (rng3.get_unique_num() != rng.get_unique_num())
    ITFAILS;
  if (rng3.get_unique_num() != rng2.get_unique_num())
    ITFAILS;
  if (!soft_equiv(w, y))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("test_stream passed");
}

//---------------------------------------------------------------------------//
void test_alias(UnitTest &ut) {
  // Create four Counter_RNG2s; rng and rng2 are identical, and rng, rng2, and
  // rng3 have the same stream number.
  uint64_t const streamnum = 0x20202020;
  uint64_t const time_step = 799;
  Counter_RNG2 rng(0x1111, streamnum, time_step);
  Counter_RNG2 rng2(0x1111, streamnum, time_step);
  Counter_RNG2 rng3(0x2222, streamnum, time_step);
  Counter_RNG2 rng4(0x3333, streamnum+1, time_step);

  if (rng.get_stream_number() != rng2.get_stream_number())
    ITFAILS;
  if (rng.get_stream_number() != rng3.get_stream_number())
    ITFAILS;
  if (rng.get_stream_number() == rng4.get_stream_number())
    ITFAILS;
  if (rng2.get_stream_number() != rng3.get_stream_number())
    ITFAILS;
  if (rng2.get_stream_number() == rng4.get_stream_number())
    ITFAILS;
  if (rng3.get_stream_number() == rng4.get_stream_number())
    ITFAILS;
  if (rng.get_unique_num() != rng2.get_unique_num())
    ITFAILS;
  if (rng.get_unique_num() == rng3.get_unique_num())
    ITFAILS;
  if (rng.get_unique_num() == rng4.get_unique_num())
    ITFAILS;
  if (rng2.get_unique_num() == rng3.get_unique_num())
    ITFAILS;
  if (rng2.get_unique_num() == rng4.get_unique_num())
    ITFAILS;
  if (rng3.get_unique_num() == rng4.get_unique_num())
    ITFAILS;
  if (rng != rng2)
    ITFAILS;
  if (rng == rng3)
    ITFAILS;
  if (rng == rng4)
    ITFAILS;
  if (rng2 == rng3)
    ITFAILS;
  if (rng2 == rng4)
    ITFAILS;
  if (rng3 == rng4)
    ITFAILS;

  // Create a Counter_RNG2_Ref from rng.
  Counter_RNG2_Ref ref(rng.ref());

  if (ref.get_stream_number() != rng.get_stream_number())
    ITFAILS;
  if (ref.get_stream_number() != rng2.get_stream_number())
    ITFAILS;
  if (ref.get_stream_number() != rng3.get_stream_number())
    ITFAILS;
  if (ref.get_stream_number() == rng4.get_stream_number())
    ITFAILS;
  if (ref.get_unique_num() != rng.get_unique_num())
    ITFAILS;
  if (ref.get_unique_num() != rng2.get_unique_num())
    ITFAILS;
  if (ref.get_unique_num() == rng3.get_unique_num())
    ITFAILS;
  if (ref.get_unique_num() == rng4.get_unique_num())
    ITFAILS;
  if (!ref.is_alias_for(rng))
    ITFAILS;
  if (ref.is_alias_for(rng2))
    ITFAILS;
  if (ref.is_alias_for(rng3))
    ITFAILS;
  if (ref.is_alias_for(rng4))
    ITFAILS;

  // Generate a random double (and advance the stream) from rng via ref.
  double x = ref.ran();

  if (ref.get_stream_number() != rng.get_stream_number())
    ITFAILS;
  if (ref.get_stream_number() != rng2.get_stream_number())
    ITFAILS;
  if (ref.get_stream_number() != rng3.get_stream_number())
    ITFAILS;
  if (ref.get_stream_number() == rng4.get_stream_number())
    ITFAILS;
  if (ref.get_unique_num() != rng.get_unique_num())
    ITFAILS;
  if (ref.get_unique_num() != rng2.get_unique_num())
    ITFAILS;
  if (ref.get_unique_num() == rng3.get_unique_num())
    ITFAILS;
  if (ref.get_unique_num() == rng4.get_unique_num())
    ITFAILS;
  if (!ref.is_alias_for(rng))
    ITFAILS;
  if (ref.is_alias_for(rng2))
    ITFAILS;
  if (ref.is_alias_for(rng3))
    ITFAILS;
  if (ref.is_alias_for(rng4))
    ITFAILS;

  // Invoking ref.ran should have altered rng; it should still have the same
  // stream number as rng2 and rng3, but it should be identical to none of them.
  if (rng.get_stream_number() != rng2.get_stream_number())
    ITFAILS;
  if (rng.get_stream_number() != rng3.get_stream_number())
    ITFAILS;
  if (rng.get_stream_number() == rng4.get_stream_number())
    ITFAILS;
  if (rng.get_unique_num() != rng2.get_unique_num())
    ITFAILS;
  if (rng.get_unique_num() == rng3.get_unique_num())
    ITFAILS;
  if (rng.get_unique_num() == rng4.get_unique_num())
    ITFAILS;
  if (rng == rng2)
    ITFAILS;
  if (rng == rng3)
    ITFAILS;
  if (rng == rng4)
    ITFAILS;

  // Create a bare data array that should match rng.
  vector<uint64_t> data(CBRNG_DATA_SIZE);
  data[0] = 1;
  data[1] = streamnum;
  data[2] = static_cast<uint64_t>(0x1111);
  data[3] = time_step;

  if (!std::equal(rng.begin(), rng.end(), data.begin()))
    ITFAILS;

  // Create a Counter_RNG2_Ref from a bare data array.
  data[0] = 0;
  data[1] = streamnum;
  data[2] = static_cast<uint64_t>(0x2222);
  data[3] = time_step;
  Counter_RNG2_Ref ref2(&data[0], &data[0] + CBRNG_DATA_SIZE);

  // ref2 should have the same stream number as rng, rng2, and rng3 but
  // shouldn't be an alias for any of them.
  if (ref2.get_stream_number() != rng.get_stream_number())
    ITFAILS;
  if (ref2.get_stream_number() != rng2.get_stream_number())
    ITFAILS;
  if (ref2.get_stream_number() != rng3.get_stream_number())
    ITFAILS;
  if (ref2.get_stream_number() == rng4.get_stream_number())
    ITFAILS;
  if (ref2.get_unique_num() == rng.get_unique_num())
    ITFAILS;
  if (ref2.get_unique_num() == rng2.get_unique_num())
    ITFAILS;
  if (ref2.get_unique_num() != rng3.get_unique_num())
    ITFAILS;
  if (ref2.get_unique_num() == rng4.get_unique_num())
    ITFAILS;
  if (ref2.is_alias_for(rng))
    ITFAILS;
  if (ref2.is_alias_for(rng2))
    ITFAILS;
  if (ref2.is_alias_for(rng3))
    ITFAILS;
  if (ref2.is_alias_for(rng4))
    ITFAILS;

  // Generate a random double from ref.
  double y = ref2.ran();

  // The underlying data array should have changed.
  if (data[0] != 1)
    ITFAILS;
  if (data[2] != static_cast<uint64_t>(0x2222))
    ITFAILS;
  if (data[1] != 0x20202020)
    ITFAILS;
  if (data[3] != time_step)
    ITFAILS;
  if (soft_equiv(y, x))
    ITFAILS;

  // Generate a random double from rng3; it should match the one from ref2.
  double z = rng3.ran();

  if (!soft_equiv(z, y))
    ITFAILS;
  if (soft_equiv(z, x))
    ITFAILS;

// Try to create a Counter_RNG2_Ref with a data array that's too short.
// 1. Only test exceptions if DbC is enabled.
// 2. However, do not run these tests if no-throw DbC is enabled (DBC & 8)
#ifdef REQUIRE_ON
#if !(DBC & 8)
  bool caught = false;
  try {
    Counter_RNG2_Ref ref3(&data[0], &data[0] + CBRNG_DATA_SIZE - 1);
  } catch (rtt_dsxx::assertion &err) {
    cout << "Good, caught assertion: " << err.what() << endl;
    caught = true;
  }
  if (!caught)
    ITFAILS;
#endif
#endif

  if (ut.numFails == 0)
    PASSMSG("test_alias passed");
}

//---------------------------------------------------------------------------//
void test_rollover(UnitTest &ut) {
  // Create a Counter_RNG2 with a large counter value.
  // data[0] = fffffffd;
  // data[1] = 1;
  // data[2] = 0xabcd;
  // data[3] = 0xef00;
  Counter_RNG2 rng;
  uint32_t const lotta_steps(0xFFFFFFFD);
  uint32_t const seed(0xCEED);
  uint64_t const time_step(99723);
  uint64_t const stream_number(0xC0FFEE);

  rng.set_step_counter(lotta_steps);
  rng.set_seed(seed);
  rng.set_time_step(time_step);
  rng.set_stream_number(stream_number);
  // convenient lambda for testing
  auto rest_unchanged = [=](Counter_RNG2 const &rng){
    return seed == rng.get_seed() &&
           time_step == rng.get_time_step() &&
           stream_number == rng.get_stream_number();
  };

  // Increment data[0], generate a random double, and compare.
  double x = rng.ran();
  FAIL_IF_NOT(lotta_steps + 1 == rng.get_step_counter());
  FAIL_IF_NOT(rest_unchanged(rng));

  // ... and again.
  double y = rng.ran();
  if (soft_equiv(y, x))
    ITFAILS;
  FAIL_IF_NOT(lotta_steps + 2 == rng.get_step_counter());
  FAIL_IF_NOT(rest_unchanged(rng));

  // Generate another random double and verify that the counter has incremented
  // correctly.
  FAIL_IF(rng.passed_max());
  double z = rng.ran();
  if (soft_equiv(z, x))
    ITFAILS;
  if (soft_equiv(z, y))
    ITFAILS;
  FAIL_IF_NOT(0 == rng.get_step_counter());
  FAIL_IF_NOT(rng.passed_max());
  FAIL_IF_NOT(rest_unchanged(rng));

  // Repeat the test with a Counter_RNG2_Ref.
  vector<uint64_t> data(CBRNG_DATA_SIZE);
  data[0] = lotta_steps + 1;
  data[1] = stream_number;
  data[2] = seed;
  data[3] = time_step;
  Counter_RNG2_Ref ref(&data[0], &data[0] + CBRNG_DATA_SIZE);

  double y2 = ref.ran();
  if (!soft_equiv(y2, y))
    ITFAILS;
  if (data[0] != lotta_steps + 2)
    ITFAILS;
  if (data[1] != stream_number)
    ITFAILS;
  if (data[2] != seed)
    ITFAILS;
  if (data[3] != time_step)
    ITFAILS;

  double z2 = ref.ran();
  if (!soft_equiv(z2, z))
    ITFAILS;
  if (data[0] != 0)
    ITFAILS;
  if (data[1] != stream_number)
    ITFAILS;
  if (data[2] != seed)
    ITFAILS;
  if (data[3] != time_step)
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("test_rollover passed");
} // test_rollover(UnitTest &ut)

//---------------------------------------------------------------------------//
void test_spawn(UnitTest &ut) {
  // Create a generator.
  uint32_t const seed = 0xabcdef;
  uint64_t const streamnum = 0;
  uint64_t const time_step = 10001;
  Counter_RNG2 rng(seed, streamnum, time_step);

  // Compare two rngs, making sure what should be different is and what
  // shouldn't be different isn't.
  auto close_but_not_too_close = [&](Counter_RNG2 const &child,
                                     Counter_RNG2 const &parent) {
    FAIL_IF(child == parent);
    FAIL_IF(child.get_stream_number() != parent.get_stream_number());
    FAIL_IF(child.get_time_step() != parent.get_time_step());
    FAIL_IF(child.get_seed() != parent.get_seed());
    FAIL_IF(child.get_spawn_id() == parent.get_spawn_id());
    FAIL_IF(child.get_spawn_id() == 0);
    return;
  };

  // Can spawn?
  FAIL_IF_NOT(rng.can_spawn());

  // Spawn a new generator.
  Counter_RNG2 rng_child1;
  rng.spawn(rng_child1, 1ull);

  // child should have same stream number, time step, and seed. It should
  // should have different spawn id. Its step counter should be zero.
  close_but_not_too_close(rng_child1, rng);

  // Create a reference to rng, and spawn from the reference.
  Counter_RNG2_Ref rng_ref(rng.ref());
  if (!rng_ref.is_alias_for(rng))
    ITFAILS;

  FAIL_IF_NOT(rng_ref.can_spawn());
  Counter_RNG2 rng_child2;
  rng_ref.spawn(rng_child2,2ull);

  // rng_child2 should have the same stream number as rng and rng_child1 but
  // should not be identical to either previous generator.
  close_but_not_too_close(rng_child2, rng);
  FAIL_IF(rng_child2 == rng_child1);

  // Spawn a generator from rng_child1.
  FAIL_IF_NOT(rng_child1.can_spawn());
  Counter_RNG2 rng_grandchild1;
  rng_child1.spawn(rng_grandchild1,3);
  close_but_not_too_close(rng_grandchild1, rng_child1);
  close_but_not_too_close(rng_grandchild1, rng);

  // Spawn a generator from rng_child2.
  FAIL_IF_NOT(rng_child2.can_spawn());
  Counter_RNG2 rng_grandchild2;
  rng_child2.spawn(rng_grandchild2,4);
  close_but_not_too_close(rng_grandchild2, rng_child2);

  // ensure that RNG::can_spawn returns false after the right number of gens
  uint32_t const mx_gens = max_gens();
  std::vector<Counter_RNG2> rngs(mx_gens+1);
  for(size_t i = 0; i < mx_gens; ++i){
    Counter_RNG2 &g = rngs[i];
    Counter_RNG2 &parent = i == 0 ? rng : rngs[i-1];

    parent.spawn(g,i+1);
    if(i < mx_gens-1){
      FAIL_IF_NOT(g.can_spawn());
    }
    else{
      FAIL_IF(g.can_spawn());
      /* should never get here--if it does, interrogate a bit */
      if(g.can_spawn()){
        printf("%s:%i i = %lu, mx_gens = %u\n",__FUNCTION__,__LINE__,i,mx_gens);
      }
    }
  }
printf("%s:%i \n",__FUNCTION__,__LINE__);
  if (ut.numFails == 0)
    PASSMSG("test_spawn passed");
}

//---------------------------------------------------------------------------//
void test_unique(UnitTest &ut) {
  // Create three identical generators.
  uint32_t const seed = 332211;
  uint64_t const streamnum = 2468;
  uint64_t const time_step = 147963;
  Counter_RNG2 rng(seed, streamnum, time_step);
  Counter_RNG2 rng2(seed, streamnum, time_step);
  Counter_RNG2 rng3(seed, streamnum, time_step);

  Counter_RNG2_Ref rng_ref(rng.ref());
  Counter_RNG2_Ref rng2_ref(rng2.ref());
  Counter_RNG2_Ref rng3_ref(rng3.ref());

  if (rng != rng2)
    ITFAILS;
  if (rng != rng3)
    ITFAILS;
  if (rng.get_unique_num() != rng2.get_unique_num())
    ITFAILS;
  if (rng.get_unique_num() != rng3.get_unique_num())
    ITFAILS;

  if (!rng_ref.is_alias_for(rng))
    ITFAILS;
  if (!rng2_ref.is_alias_for(rng2))
    ITFAILS;
  if (!rng3_ref.is_alias_for(rng3))
    ITFAILS;

  // Generate some random numbers from rng2.  The stream number and unique
  // number of rng2 should remain the same during this process.
  set<uint64_t> ids;
  ids.insert(rng.get_unique_num());

  for (int i = 0; i < 1000000; ++i) {
    rng2.ran();

    if (rng2.get_stream_number() != rng.get_stream_number())
      ITFAILS;
    if (rng2_ref.get_unique_num() != rng2.get_unique_num())
      ITFAILS;
    if (ids.find(rng2.get_unique_num()) == ids.end())
      ITFAILS;
  }

  if (ut.numFails == 0)
    PASSMSG("test_unique passed");
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    test_set_get_seed(ut);
    test_set_get_step_counter(ut);
    test_set_get_time_step(ut);
    test_set_get_stream_number(ut);
    test_set_get_spawn_id(ut);
    test_increment_CBRNG2_step_counter(ut);
    test_equality(ut);
    test_stream(ut);
    test_alias(ut);
    test_rollover(ut);
    test_spawn(ut);
    test_unique(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstCounter_RNG2.cc
//---------------------------------------------------------------------------//
