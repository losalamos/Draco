//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   rng/test/tstCounter_RNG.cc
 * \author Peter Ahrens
 * \date   Fri Aug 3 16:53:23 2012
 * \brief  Counter_RNG tests.
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "rng/Counter_RNG.hh"
#include <set>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_rng;

//----------------------------------------------------------------------------//
// TESTS
//----------------------------------------------------------------------------//

void test_equality(UnitTest &ut) {
  // Create a Counter_RNG by specifying a seed and stream number.
  uint32_t seed = 1;
  uint64_t streamnum = 2;
  Counter_RNG rng(seed, streamnum);

  if (rng.get_num() != streamnum)
    ITFAILS;
  if (rng.size() != CBRNG_DATA_SIZE)
    ITFAILS;
  if (rng.size_bytes() != CBRNG_DATA_SIZE * sizeof(uint64_t))
    ITFAILS;
  if (rng != rng)
    ITFAILS;

  // Create another Counter_RNG with a different seed.
  seed = 2;
  Counter_RNG rng2(seed, streamnum);

  // rng2's stream number should match rng's, but the two generators should not
  // be identical.
  if (rng2.get_num() != streamnum)
    ITFAILS;
  if (rng2.get_num() != rng.get_num())
    ITFAILS;
  if (rng2.get_unique_num() == rng.get_unique_num())
    ITFAILS;
  if (rng2 == rng)
    ITFAILS;
  if (rng2 != rng2)
    ITFAILS;

  // Create another Counter_RNG with a different stream number.
  seed = 1;
  streamnum = 3;
  Counter_RNG rng3(seed, streamnum);

  // rng3 should be different from the previous two generators.
  if (rng3.get_num() != streamnum)
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

  // Create another Counter_RNG with the original seed and stream number.
  streamnum = 2;
  Counter_RNG rng4(seed, streamnum);

  // rng4 should be equal to rng but different from rng2 and rng3.
  if (rng4.get_num() != streamnum)
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

  // Create a Counter_RNG from a data array.
  vector<uint64_t> data(CBRNG_DATA_SIZE);
  data[0] = 1234;
  data[1] = 5678;
  data[2] = 9012;
  data[3] = 3456;
  Counter_RNG rng5(&data[0], &data[0] + CBRNG_DATA_SIZE);

  streamnum = data[2];
  if (rng5.get_num() != streamnum)
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

  // Create a Counter_RNG from a data array that should match rng and rng4.
  data[0] = 0;
  data[1] = static_cast<uint64_t>(1) << 32;
  data[2] = 2;
  data[3] = 0;
  Counter_RNG rng6(&data[0], &data[0] + CBRNG_DATA_SIZE);

  streamnum = data[2];
  if (rng6.get_num() != streamnum)
    ITFAILS;
  if (rng6.get_unique_num() != rng.get_unique_num())
    ITFAILS;
  if (rng6.get_unique_num() == rng2.get_unique_num())
    ITFAILS;
  if (rng6.get_unique_num() == rng3.get_unique_num())
    ITFAILS;
  if (rng6.get_unique_num() != rng4.get_unique_num())
    ITFAILS;
  if (rng6.get_unique_num() == rng5.get_unique_num())
    ITFAILS;
  if (rng6 != rng)
    ITFAILS;
  if (rng6 == rng2)
    ITFAILS;
  if (rng6 == rng3)
    ITFAILS;
  if (rng6 != rng4)
    ITFAILS;
  if (rng6 == rng5)
    ITFAILS;
  if (rng6 != rng6)
    ITFAILS;

// Try to create a Counter_RNG from a data array that's too short.
// 1. Only test exceptions if DbC is enabled.
// 2. However, do not run these tests if no-throw DbC is enabled (DBC & 8)
#ifdef REQUIRE_ON
#if !(DBC & 8)
  bool caught = false;
  try {
    Counter_RNG rng7(&data[0], &data[0] + CBRNG_DATA_SIZE - 1);
  } catch (rtt_dsxx::assertion &err) {
    cout << "Good, caught assertion: " << err.what() << endl;
    caught = true;
  }
  if (!caught)
    ITFAILS;
#endif
#endif

  // Test for equality using iterators.
  if (!std::equal(rng6.begin(), rng6.end(), rng.begin()))
    ITFAILS;
  if (std::equal(rng6.begin(), rng6.end(), rng2.begin()))
    ITFAILS;
  if (std::equal(rng6.begin(), rng6.end(), rng3.begin()))
    ITFAILS;
  if (!std::equal(rng6.begin(), rng6.end(), rng4.begin()))
    ITFAILS;
  if (std::equal(rng6.begin(), rng6.end(), rng5.begin()))
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("test_equality passed");
}

//----------------------------------------------------------------------------//
void test_stream(UnitTest &ut) {
  // Create two identical Counter_RNGs.
  uint32_t seed = 0x12121212;
  uint64_t streamnum = 1234;
  Counter_RNG rng(seed, streamnum);
  Counter_RNG rng2(seed, streamnum);

  if (rng != rng2)
    ITFAILS;

  // Generate a random double (and advance the stream) from rng.
  double x = rng.ran();

  // rng and rng2 should no longer match, but their stream numbers and unique
  // identifiers should be the same.
  if (rng == rng2)
    ITFAILS;
  if (rng.get_num() != streamnum)
    ITFAILS;
  if (rng.get_num() != rng2.get_num())
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
  if (rng.get_num() != streamnum)
    ITFAILS;
  if (rng.get_num() != rng2.get_num())
    ITFAILS;
  if (rng.get_unique_num() != rng2.get_unique_num())
    ITFAILS;
  if (soft_equiv(x, z))
    ITFAILS;

  // Create a Counter_RNG from a data array.
  vector<uint64_t> data(CBRNG_DATA_SIZE);
  data[0] = 0;
  data[1] = static_cast<uint64_t>(seed) << 32;
  data[2] = streamnum;
  data[3] = 0;
  Counter_RNG rng3(&data[0], &data[0] + CBRNG_DATA_SIZE);

  // Initially, rng3 should exactly match neither rng nor rng2, but all three
  // should have the same stream number and "unique" identifier.
  if (!std::equal(rng3.begin(), rng3.end(), data.begin()))
    ITFAILS;
  if (rng3 == rng)
    ITFAILS;
  if (rng3 == rng2)
    ITFAILS;
  if (rng3.get_num() != streamnum)
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
  if (rng3.get_num() != streamnum)
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

//----------------------------------------------------------------------------//
void test_alias(UnitTest &ut) {
  // Create four Counter_RNGs; rng and rng2 are identical, and rng, rng2, and
  // rng3 have the same stream number.
  uint64_t streamnum = 0x20202020;
  Counter_RNG rng(0x1111, streamnum);
  Counter_RNG rng2(0x1111, streamnum);
  Counter_RNG rng3(0x2222, streamnum);
  ++streamnum;
  Counter_RNG rng4(0x3333, streamnum);

  if (rng.get_num() != rng2.get_num())
    ITFAILS;
  if (rng.get_num() != rng3.get_num())
    ITFAILS;
  if (rng.get_num() == rng4.get_num())
    ITFAILS;
  if (rng2.get_num() != rng3.get_num())
    ITFAILS;
  if (rng2.get_num() == rng4.get_num())
    ITFAILS;
  if (rng3.get_num() == rng4.get_num())
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

  // Create a Counter_RNG_Ref from rng.
  Counter_RNG_Ref ref(rng.ref());

  if (ref.get_num() != rng.get_num())
    ITFAILS;
  if (ref.get_num() != rng2.get_num())
    ITFAILS;
  if (ref.get_num() != rng3.get_num())
    ITFAILS;
  if (ref.get_num() == rng4.get_num())
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

  if (ref.get_num() != rng.get_num())
    ITFAILS;
  if (ref.get_num() != rng2.get_num())
    ITFAILS;
  if (ref.get_num() != rng3.get_num())
    ITFAILS;
  if (ref.get_num() == rng4.get_num())
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
  if (rng.get_num() != rng2.get_num())
    ITFAILS;
  if (rng.get_num() != rng3.get_num())
    ITFAILS;
  if (rng.get_num() == rng4.get_num())
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
  data[1] = static_cast<uint64_t>(0x1111) << 32;
  data[2] = 0x20202020;
  data[3] = 0;

  if (!std::equal(rng.begin(), rng.end(), data.begin()))
    ITFAILS;

  // Create a Counter_RNG_Ref from a bare data array.
  data[0] = 0;
  data[1] = static_cast<uint64_t>(0x2222) << 32;
  data[2] = 0x20202020;
  data[3] = 0;
  Counter_RNG_Ref ref2(&data[0], &data[0] + CBRNG_DATA_SIZE);

  // ref2 should have the same stream number as rng, rng2, and rng3 but
  // shouldn't be an alias for any of them.
  if (ref2.get_num() != rng.get_num())
    ITFAILS;
  if (ref2.get_num() != rng2.get_num())
    ITFAILS;
  if (ref2.get_num() != rng3.get_num())
    ITFAILS;
  if (ref2.get_num() == rng4.get_num())
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
  if (data[1] != static_cast<uint64_t>(0x2222) << 32)
    ITFAILS;
  if (data[2] != 0x20202020)
    ITFAILS;
  if (data[3] != 0)
    ITFAILS;
  if (soft_equiv(y, x))
    ITFAILS;

  // Generate a random double from rng3; it should match the one from ref2.
  double z = rng3.ran();

  if (!soft_equiv(z, y))
    ITFAILS;
  if (soft_equiv(z, x))
    ITFAILS;

// Try to create a Counter_RNG_Ref with a data array that's too short.
// 1. Only test exceptions if DbC is enabled.
// 2. However, do not run these tests if no-throw DbC is enabled (DBC & 8)
#ifdef REQUIRE_ON
#if !(DBC & 8)
  bool caught = false;
  try {
    Counter_RNG_Ref ref3(&data[0], &data[0] + CBRNG_DATA_SIZE - 1);
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

//----------------------------------------------------------------------------//
void test_rollover(UnitTest &ut) {
  // Create a Counter_RNG with a large counter value.
  vector<uint64_t> data(CBRNG_DATA_SIZE);
  data[0] = 0xfffffffffffffffd;
  data[1] = 1;
  data[2] = 0xabcd;
  data[3] = 0xef00;
  Counter_RNG rng(&data[0], &data[0] + CBRNG_DATA_SIZE);

  // Increment data[0], generate a random double, and compare.
  ++data[0];
  double x = rng.ran();
  if (!std::equal(rng.begin(), rng.end(), data.begin()))
    ITFAILS;

  // ... and again.
  ++data[0];
  double y = rng.ran();
  if (soft_equiv(y, x))
    ITFAILS;
  if (!std::equal(rng.begin(), rng.end(), data.begin()))
    ITFAILS;

  // Generate another random double and verify that the counter has incremented
  // correctly.
  data[0] = 0;
  data[1] = 2;
  double z = rng.ran();
  if (soft_equiv(z, x))
    ITFAILS;
  if (soft_equiv(z, y))
    ITFAILS;
  if (!std::equal(rng.begin(), rng.end(), data.begin()))
    ITFAILS;

  // Repeat the test with a Counter_RNG_Ref.
  data[0] = 0xfffffffffffffffe;
  data[1] = 1;
  Counter_RNG_Ref ref(&data[0], &data[0] + CBRNG_DATA_SIZE);

  double y2 = ref.ran();
  if (!soft_equiv(y2, y))
    ITFAILS;
  if (data[0] != 0xffffffffffffffff)
    ITFAILS;
  if (data[1] != 1)
    ITFAILS;
  if (data[2] != 0xabcd)
    ITFAILS;
  if (data[3] != 0xef00)
    ITFAILS;

  double z2 = ref.ran();
  if (!soft_equiv(z2, z))
    ITFAILS;
  if (data[0] != 0)
    ITFAILS;
  if (data[1] != 2)
    ITFAILS;
  if (data[2] != 0xabcd)
    ITFAILS;
  if (data[3] != 0xef00)
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("test_rollover passed");
}

//----------------------------------------------------------------------------//
void test_spawn(UnitTest &ut) {
  // Create a generator.
  uint32_t seed = 0xabcdef;
  uint64_t streamnum = 0;
  Counter_RNG rng(seed, streamnum);

  // Spawn a new generator.
  Counter_RNG rng_child1;
  rng.spawn(rng_child1);

  // rng_child1 should have the same stream number as rng but should not be
  // identical to rng.
  if (rng_child1.get_num() != streamnum)
    ITFAILS;
  if (rng.get_num() != streamnum)
    ITFAILS;
  if (rng_child1 == rng)
    ITFAILS;

  // Create a reference to rng, and spawn from the reference.
  Counter_RNG_Ref rng_ref(rng.ref());
  if (!rng_ref.is_alias_for(rng))
    ITFAILS;

  Counter_RNG rng_child2;
  rng_ref.spawn(rng_child2);

  // rng_child2 should have the same stream number as rng and rng_child1 but
  // should not be identical to either previous generator.
  if (rng_child2.get_num() != streamnum)
    ITFAILS;
  if (rng_child1.get_num() != streamnum)
    ITFAILS;
  if (rng.get_num() != streamnum)
    ITFAILS;
  if (rng_child2 == rng_child1)
    ITFAILS;
  if (rng_child2 == rng)
    ITFAILS;
  if (rng_child1 == rng)
    ITFAILS;

  // Spawn a generator from rng_child1.
  Counter_RNG rng_grandchild1;
  rng_child1.spawn(rng_grandchild1);

  if (rng_grandchild1.get_num() != streamnum)
    ITFAILS;
  if (rng_child2.get_num() != streamnum)
    ITFAILS;
  if (rng_child1.get_num() != streamnum)
    ITFAILS;
  if (rng.get_num() != streamnum)
    ITFAILS;
  if (rng_grandchild1 == rng_child2)
    ITFAILS;
  if (rng_grandchild1 == rng_child1)
    ITFAILS;
  if (rng_grandchild1 == rng)
    ITFAILS;
  if (rng_child2 == rng_child1)
    ITFAILS;
  if (rng_child2 == rng)
    ITFAILS;
  if (rng_child1 == rng)
    ITFAILS;

  // Spawn a generator from rng_child2.
  Counter_RNG rng_grandchild2;
  rng_child2.spawn(rng_grandchild2);

  if (rng_grandchild2.get_num() != streamnum)
    ITFAILS;
  if (rng_grandchild1.get_num() != streamnum)
    ITFAILS;
  if (rng_child2.get_num() != streamnum)
    ITFAILS;
  if (rng_child1.get_num() != streamnum)
    ITFAILS;
  if (rng.get_num() != streamnum)
    ITFAILS;
  if (rng_grandchild2 == rng_grandchild1)
    ITFAILS;
  if (rng_grandchild2 == rng_child2)
    ITFAILS;
  if (rng_grandchild2 == rng_child1)
    ITFAILS;
  if (rng_grandchild2 == rng)
    ITFAILS;
  if (rng_grandchild1 == rng_child2)
    ITFAILS;
  if (rng_grandchild1 == rng_child1)
    ITFAILS;
  if (rng_grandchild1 == rng)
    ITFAILS;
  if (rng_child2 == rng_child1)
    ITFAILS;
  if (rng_child2 == rng)
    ITFAILS;
  if (rng_child1 == rng)
    ITFAILS;

  // Create two identical generators.
  Counter_RNG original(seed, streamnum);
  Counter_RNG parent(seed, streamnum);

  // Repeatedly spawn from parent.  See how long it takes to create a duplicate
  // of a previous generator.
  set<uint64_t> spawn_id;
  spawn_id.insert(original.begin()[3]);

  unsigned int expected_period = 0;
  for (unsigned int i = 0; i < 8 * sizeof(uint64_t); ++i)
    expected_period += i;

  unsigned int timeout = 16384;
  unsigned int gen = 0;
  for (gen = 0; gen < timeout; ++gen) {
    Counter_RNG child;
    parent.spawn(child);

    // The child generator should always have the same stream number as the
    // parent and the original but should never be identical to either.
    if (child.get_num() != streamnum)
      ITFAILS;
    if (parent.get_num() != streamnum)
      ITFAILS;
    if (child == parent)
      ITFAILS;
    if (child == original)
      ITFAILS;
    if (parent == original)
      ITFAILS;

    // Look for both parent and child identifiers in the spawn_id set.
    uint64_t parent_id = parent.begin()[3];
    uint64_t child_id = child.begin()[3];

    if (spawn_id.find(parent_id) != spawn_id.end()) {
      std::cout << "Found parent among previous generators after " << gen
                << " generations." << std::endl;
      break;
    }

    if (spawn_id.find(child_id) != spawn_id.end()) {
      std::cout << "Found child among previous generators after " << gen
                << " generations." << std::endl;
      break;
    }

    // Insert the parent and child identifiers and continue.
    spawn_id.insert(parent_id);
    spawn_id.insert(child_id);
  }
  if (gen != expected_period)
    ITFAILS;

  // Go again from that parent, this time through a reference.  How long until
  // it repeats this time?
  Counter_RNG_Ref parent_ref(parent.ref());
  if (!parent_ref.is_alias_for(parent))
    ITFAILS;

  spawn_id.clear();
  spawn_id.insert(original.begin()[3]);
  spawn_id.insert(parent.begin()[3]);

  // Start from generation 1, to include the generational difference in this
  // experiment between original and parent.
  for (gen = 1; gen < timeout; ++gen) {
    Counter_RNG child;
    parent_ref.spawn(child);
    if (!parent_ref.is_alias_for(parent))
      ITFAILS;

    // The child generator should always have the same stream number as the
    // parent and the original but should never be identical to either.
    if (child.get_num() != streamnum)
      ITFAILS;
    if (parent.get_num() != streamnum)
      ITFAILS;
    if (child == parent)
      ITFAILS;
    if (child == original)
      ITFAILS;
    if (parent == original)
      ITFAILS;

    // Look for both parent and child identifiers in the spawn_id set.
    uint64_t parent_id = parent.begin()[3];
    uint64_t child_id = child.begin()[3];

    if (spawn_id.find(parent_id) != spawn_id.end()) {
      std::cout << "Found parent among previous generators after " << gen
                << " generations." << std::endl;
      break;
    }

    if (spawn_id.find(child_id) != spawn_id.end()) {
      std::cout << "Found child among previous generators after " << gen
                << " generations." << std::endl;
      break;
    }

    // Insert the parent and child identifiers and continue.
    spawn_id.insert(parent_id);
    spawn_id.insert(child_id);
  }
  if (gen != expected_period)
    ITFAILS;

  // Repeat the experiment, but this time use the first child as the starting
  // parent.
  Counter_RNG child;
  parent.spawn(child);

  if (child.get_num() != streamnum)
    ITFAILS;
  if (parent.get_num() != streamnum)
    ITFAILS;
  if (child == parent)
    ITFAILS;
  if (child == original)
    ITFAILS;
  if (parent == original)
    ITFAILS;

  spawn_id.clear();
  spawn_id.insert(original.begin()[3]);
  spawn_id.insert(parent.begin()[3]);
  spawn_id.insert(child.begin()[3]);

  expected_period = 1;
  for (unsigned int i = 0; i < 8 * sizeof(uint64_t) - 2; ++i)
    expected_period += i;

  for (gen = 0; gen < timeout; ++gen) {
    Counter_RNG grandchild;
    child.spawn(grandchild);

    // The grandchild generator should always have the same stream number as its
    // parent, grandparent, and the original generator but should never be
    // identical to any of them.
    if (grandchild.get_num() != streamnum)
      ITFAILS;
    if (child.get_num() != streamnum)
      ITFAILS;
    if (parent.get_num() != streamnum)
      ITFAILS;
    if (grandchild == child)
      ITFAILS;
    if (grandchild == parent)
      ITFAILS;
    if (grandchild == original)
      ITFAILS;

    // Look for both child and grandchild identifiers in the spawn_id set.
    uint64_t child_id = child.begin()[3];
    uint64_t grandchild_id = grandchild.begin()[3];

    if (spawn_id.find(child_id) != spawn_id.end()) {
      std::cout << "Found child among previous generators after " << gen
                << " generations." << std::endl;
      break;
    }

    if (spawn_id.find(grandchild_id) != spawn_id.end()) {
      std::cout << "Found grandchild among previous generators after " << gen
                << " generations." << std::endl;
      break;
    }

    // Insert the child and grandchild identifiers and continue.
    spawn_id.insert(child_id);
    spawn_id.insert(grandchild_id);
  }
  if (gen != expected_period)
    ITFAILS;

  if (ut.numFails == 0)
    PASSMSG("test_spawn passed");
}

//----------------------------------------------------------------------------//
void test_unique(UnitTest &ut) {
  // Create three identical generators.
  uint32_t seed = 332211;
  uint64_t streamnum = 2468;
  Counter_RNG rng(seed, streamnum);
  Counter_RNG rng2(seed, streamnum);
  Counter_RNG rng3(seed, streamnum);

  Counter_RNG_Ref rng_ref(rng.ref());
  Counter_RNG_Ref rng2_ref(rng2.ref());
  Counter_RNG_Ref rng3_ref(rng3.ref());

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

    if (rng2.get_num() != rng.get_num())
      ITFAILS;
    if (rng2_ref.get_unique_num() != rng2.get_unique_num())
      ITFAILS;
    if (ids.find(rng2.get_unique_num()) == ids.end())
      ITFAILS;
  }

  // Spawn some generators from rng3.  While all spawned generators should have
  // the same stream number as the original parent, their unique numbers should
  // differ.
  unsigned int expected_period = 0;
  for (unsigned int i = 0; i < 8 * sizeof(uint64_t); ++i)
    expected_period += i;

  for (unsigned int i = 0; i < expected_period; ++i) {
    Counter_RNG rng3_child;
    rng3.spawn(rng3_child);

    if (rng3.get_num() != rng.get_num())
      ITFAILS;
    if (rng3_child.get_num() != rng.get_num())
      ITFAILS;

    if (rng_ref.get_unique_num() != rng.get_unique_num())
      ITFAILS;
    if (rng3_ref.get_unique_num() != rng3.get_unique_num())
      ITFAILS;

    if (ids.find(rng3.get_unique_num()) != ids.end())
      ITFAILS;
    if (ids.find(rng3_child.get_unique_num()) != ids.end())
      ITFAILS;

    // Insert new unique identifiers and continue.
    ids.insert(rng3.get_unique_num());
    ids.insert(rng3_child.get_unique_num());
  }

  if (ut.numFails == 0)
    PASSMSG("test_unique passed");
}

//----------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    test_equality(ut);
    test_stream(ut);
    test_alias(ut);
    test_rollover(ut);
    test_spawn(ut);
    test_unique(ut);
  }
  UT_EPILOG(ut);
}

//----------------------------------------------------------------------------//
// end of tstCounter_RNG.cc
//----------------------------------------------------------------------------//
