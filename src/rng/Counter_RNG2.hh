//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Counter_RNG2.hh
 * \brief  (Re)Declaration of class Counter_RNG2.
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved */
//---------------------------------------------------------------------------//

#ifndef Counter_RNG2_hh
#define Counter_RNG2_hh

#include "rng/config.h"

#ifdef _MSC_FULL_VER
// Engines have multiple copy constructors, quite legal C++, disable MSVC
// complaint.
#pragma warning(disable : 4521)
#endif

#if defined(__ICC)
// Suppress Intel's "unrecognized preprocessor directive" warning, triggered by
// use of #warning in Random123/features/sse.h.
#pragma warning disable 11
#endif

#if defined(__GNUC__) && !defined(__clang__)

/*
#if (RNG_GNUC_VERSION >= 40204) && !defined(__ICC) && !defined(NVCC)
// Suppress GCC's "unused parameter" warning, about lhs and rhs in sse.h, and an
// "unused local typedef" warning, from a pre-C++11 implementation of a static
// assertion in compilerfeatures.h.
*/
#pragma GCC diagnostic push
#if (RNG_GNUC_VERSION >= 70000)
#pragma GCC diagnostic ignored "-Wexpansion-to-defined"
#endif
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wexpansion-to-defined"
#endif

#include "Random123/threefry.h"
#include "uniform.hpp"

#ifdef __clang__
// Restore clang diagnostics to previous state.
#pragma clang diagnostic pop
#endif

/* #if (RNG_GNUC_VERSION >= 40600) */
#if defined(__GNUC__) && !defined(__clang__)
/* && (RNG_GNUC_VERSION >= 70000) */
// Restore GCC diagnostics to previous state.
#pragma GCC diagnostic pop
#endif

#include "ds++/Data_Table.hh"
#include <algorithm>
#include <limits>
#include <type_traits>

namespace rtt_rng {

// Forward declaration.
class Counter_RNG2;

// Select a particular counter-based random number generator from Random123.
typedef r123::Threefry2x64 CBRNG;

// Counter and key types.
typedef CBRNG::ctr_type ctr_type;
typedef CBRNG::key_type key_type;

#define CBRNG_DATA_SIZE 4

// These implementation details get moved out here because of the
// circularity of Counter_RNG2_Ref and Counter_RNG2. Once the former is
// eliminated, these masks and getters will move back into Counter_RNG2.
static const uint64_t SEED_MASK = 0x00000000FFFFFFFF;

static const uint64_t STEP_MASK = 0x00000000FFFFFFFF;

static const uint64_t SPAWN_MASK = ~STEP_MASK;

constexpr uint32_t max_gens() { return 4; }

/* These getter functions control the indexing of the RNG's data array.
 * Stateless functions permit calling them on an RNG's data, instead of the
 * RNG object itself. Why was that a good idea? Because it allowed one
 * to spawn from an RNG_Ref or from an RNG. This is a component of the
 * NR system. It should go away when NR is retired.
 */
uint32_t get_seed_impl(ctr_type::value_type const *const rng_data) {
  static_assert(std::is_same<ctr_type::value_type,uint64_t>::value,
    "Implementation relies on known integer size!");
  uint64_t const val64 = rng_data[2] & SEED_MASK;
  return static_cast<uint32_t>(val64);
}

uint32_t get_step_counter_impl(ctr_type::value_type const *const rng_data) {
  static_assert(std::is_same<ctr_type::value_type,uint64_t>::value,
    "Implementation relies on known integer size!");
  uint64_t const val64 = rng_data[0] & STEP_MASK;
  return static_cast<uint32_t>(val64);
}

uint64_t get_time_step_impl(ctr_type::value_type const *const rng_data) {
  static_assert(std::is_same<ctr_type::value_type,uint64_t>::value,
    "Implementation relies on known integer size!");
  return rng_data[3];
}

uint64_t get_stream_number_impl(ctr_type::value_type const *const rng_data) {
  static_assert(std::is_same<ctr_type::value_type,uint64_t>::value,
    "Implementation relies on known integer size!");
  return rng_data[1];
}

uint32_t get_spawn_id_impl(ctr_type::value_type const *const rng_data) {
  static_assert(std::is_same<ctr_type::value_type,uint64_t>::value,
    "Implementation relies on known integer size!");
  uint64_t const val64 = rng_data[0] & SPAWN_MASK;
  return static_cast<uint32_t>(val64 >> 32);
}

/**\brief get the RNG's spawn generation. 0 means this generator is
 * 'original' (no parents), 1 -> 1 parent gen, 2 -> 2 parent gens etc.
 */
uint32_t get_spawn_gen_impl(ctr_type::value_type const *const rng_data) {
  uint32_t const sid = get_spawn_id_impl(rng_data);
  if(0 == sid){
    return 0;
  }
  if(0xFF >= sid){
    return 1;
  }
  if(0xFFFF >= sid){
    return 2;
  }
  if(0xFFFFFF >= sid){
    return 3;
  }
  return 4;
}

/**\brier Can a generator spawn? */
bool can_spawn_impl(ctr_type::value_type const *const rng_data){
  uint32_t const spawn_gen = get_spawn_gen_impl(rng_data);
  uint32_t const max_gen = max_gens();
  return spawn_gen < max_gen;
}

namespace // anonymous
{

//---------------------------------------------------------------------------//
/*! \brief Generate a nearly-unique identifier.
 *
 * Given a pointer to RNG state data, this function generates a 64-bit
 * identifier unique to this generator but not to the specific position of its
 * RNG stream.  In other words, the identifier associated with a given generator
 * will not change as random numbers are generated from it.  However, this
 * insensitivity to the specific stream position also means that repeated
 * spawning will eventually produce two generators with the same identifier.
 *
 * This function simply applies the chosen counter-based RNG to a shuffled
 * version of the RNG seed, stream number, and spawn indicator and then returns
 * the lower 64 bits of the result.
 */
static inline uint64_t _get_unique_num(const ctr_type::value_type *const data) {
  CBRNG hash;
  uint64_t const sp64 = static_cast<uint64_t>(get_spawn_id_impl(data));
  const ctr_type ctr = {{data[1], sp64}};
  const key_type key = {{data[3], data[2]}};
  const ctr_type result = hash(ctr, key);
  return result[0];
}

//---------------------------------------------------------------------------//
/*! \brief Increment the counter, rolling over after 32b exhausted. */
inline uint64_t increment_CBRNG2_step_counter( uint64_t in){
  uint64_t spawn64 = in & (~STEP_MASK);
  uint64_t const ctr = in & STEP_MASK;
  uint64_t const new_ctr = ctr + 1;
  return spawn64 | (new_ctr & STEP_MASK);
}


//---------------------------------------------------------------------------//
/*! \brief Generate a random double.
 *
 * Given a pointer to RNG state data, this function returns a random double in
 * the open interval (0, 1)---i.e., excluding the endpoints.
 */
static inline double _ran(ctr_type::value_type *const data) {
  CBRNG rng;

  // Assemble a counter from the first two elements in data.
  ctr_type ctr = {{data[0], data[1]}};

  // Assemble a key from the last two elements in data.
  const key_type key = {{data[2], data[3]}};

  // Invoke the counter-based rng.
  const ctr_type result = rng(ctr, key);

  // Increment the counter.
  uint64_t new_ctr = increment_CBRNG2_step_counter(data[0]);
  data[0] = new_ctr;

  // Convert the first 64 bits of the RNG output into a double-precision value
  // in the open interval (0, 1) and return it.
  return r123::u01fixedpt<double, ctr_type::value_type>(result[0]);
}

} // namespace

//===========================================================================//
/*!
 * \class Counter_RNG2_Ref
 * \brief A reference to a Counter_RNG2.
 * \deprecated{The Counter_RNG2_Ref will be removed when NR is removed.}
 *
 * Counter_RNG2_Ref provides an interface to a counter-based random number
 * generator from the Random123 library from D. E. Shaw Research
 * (http://www.deshawresearch.com/resources_random123.html).  Unlike
 * Counter_RNG2, Counter_RNG2_Ref doesn't own its own RNG state (i.e., key and
 * counter); instead, it operates using a data block specified during
 * construction.
 */
//===========================================================================//
class Counter_RNG2_Ref {
public:
  //! Constructor.  db and de specify the extents of an RNG state block.
  Counter_RNG2_Ref(ctr_type::value_type *const db,
                  ctr_type::value_type *const de)
      : data(db, de) {
    Require(std::distance(db, de) * sizeof(ctr_type::value_type) ==
            sizeof(ctr_type) + sizeof(key_type));
  }

  //! Return a random double in the open interval (0, 1).
  double ran() const { return _ran(data.access()); }

  //! Spawn a new, independent generator from this reference.
  inline void spawn(Counter_RNG2 &new_gen, uint32_t const child_idx) const;

  //! get the stream number
  uint64_t get_stream_number() const {
    return get_stream_number_impl(data.access());
  }

  //! Return a unique identifier for this generator.
  uint64_t get_unique_num() const { return _get_unique_num(data.access()); }

  //! Is this Counter_RNG2_Ref a reference to rng?
  inline bool is_alias_for(Counter_RNG2 const &rng) const;

  //! Return the stream number
  uint64_t get_num() const { return get_stream_number_impl(data.access()); }

  //! Can the referred RNG spawn?
  bool can_spawn() const {
    return can_spawn_impl(data.access());
  }

private:
  mutable rtt_dsxx::Data_Table<ctr_type::value_type> data;
};

//===========================================================================//
/*!
 * \class Counter_RNG2
 * \brief A counter-based random-number generator.
 *
 * Counter_RNG2 provides an interface to a counter-based random number generator
 * from the Random123 library from D. E. Shaw Research
 * (http://www.deshawresearch.com/resources_random123.html).
 *
 * Counter_RNG2 changes how data are interpretted from the original Counter_RNG.
 * In this version, the RNG state is interpretted as five quantities: the time
 * step, the unique stream number, the spawn id, the user-provided seed, and
 * the actual step counter.
 *
 * The 128 bits of the counter are divied between a 64b
 * stream id, a 32b counter, and a 32b spawn id. The key consists of time step
 * (64b) and the random number seed (32b). The last 32b, er, reserved for future
 * use. The key is expected to be the same for all particles, meaning that a
 * thrifty programmer could save 16B/particle storage.
 *
 *  Counter |________|____|____|
 *            strm #   spwn ctr
 *
 *  Key |________|____|____|
 *       time stp| res|seed|
 *
 * Spawning is limited to four generations maximum, with each generation
 * admitting 256 children.
 *
 * Another change is the behavior on rollover. Previously, the counter had
 * 96b of counter space for RNs. Now, the counter only has 32b.
 *
 * Counter_RNG2_Ref is a friend of Counter_RNG2 because spawning a new generator
 * modifies both the parent and the child generator in ways that should not be
 * exposed through the public interface of Counter_RNG2.
 *
 * Similarly, Rnd_Control is a friend of Counter_RNG2 because initializing a
 * generator requires access to private data that should not be exposed through
 * the public interface.  Rnd_Control takes no responsibility for instantiating
 * Counter_RNG2s itself, and since copying Counter_RNG2s is disabled (via a
 * private copy constructor), an Rnd_Control must be able to initialize a
 * generator that was instantiated outside of its control.
 */
//===========================================================================//
class Counter_RNG2 {
  friend class Counter_RNG2_Ref;
  friend class Rnd_Control;

public:
  typedef ctr_type::const_iterator const_iterator;

  static constexpr uint32_t max_steps = std::numeric_limits<uint32_t>::max();

  /*! \brief Default constructor.
   *
   * This default constructor is invoked when a client wants to create a
   * Counter_RNG2 but delegate its initialization to an Rnd_Control object.
   */
  Counter_RNG2() {
    Require(sizeof(data) == sizeof(ctr_type) + sizeof(key_type));
    for(size_t i = 0; i < CBRNG_DATA_SIZE; ++i){
      data[i] = 0;
    }
  }

  //! Construct a Counter_RNG2 using a seed and stream number.
  Counter_RNG2(const uint32_t seed, const uint64_t streamnum,
               uint64_t const timestep) {
    for(size_t i = 0; i < CBRNG_DATA_SIZE; ++i){
      data[i] = 0;
    }
    initialize(seed, streamnum, timestep);
  }

  //! Create a new Counter_RNG2 from data.
  Counter_RNG2(const ctr_type::value_type *const begin,
               const ctr_type::value_type *const end) {
    Require(std::distance(begin, end) * sizeof(ctr_type::value_type) ==
            sizeof(ctr_type) + sizeof(key_type));

    std::copy(begin, end, data);
  }

  //! Return a random double in the interval (0, 1).
  double ran() const {
    if(get_step_counter() == max_steps){
      passed_max_ = true;
    }
    return _ran(data);
  }

  //! Spawn a new, independent generator from this one.
  void spawn(Counter_RNG2 &new_gen, uint32_t const child_idx) const {
    new_gen._spawn(data,child_idx);
  }

  //! Return a unique identifier for this generator.
  uint64_t get_unique_num() const { return _get_unique_num(data); }

  //! Return an iterator to the beginning of the state block.
  const_iterator begin() const { return data; }

  //! Return an iterator to the end of the state block.
  const_iterator end() const { return data + size(); }

  //! Test for equality.
  bool operator==(Counter_RNG2 const &rhs) const {
    return std::equal(begin(), end(), rhs.begin());
  }

  //! Test for inequality.
  bool operator!=(Counter_RNG2 const &rhs) const {
    return !std::equal(begin(), end(), rhs.begin());
  }

  //! Return a Counter_RNG2_Ref corresponding to this Counter_RNG2.
  Counter_RNG2_Ref ref() const { return Counter_RNG2_Ref(data, data + size()); }

  //! Return the size of this Counter_RNG2.
  size_t size() const { return size_bytes() / sizeof(ctr_type::value_type); }

  //! Return the size of this Counter_RNG2 in bytes.
  size_t size_bytes() const { return sizeof(data); }

  //! Has this counter passed the maximum number of steps?
  bool passed_max() const {
    return passed_max_;
  }

  /* Field guide to layout of bits and elements in 'data' array
   *          | data[1]| data[0] |
   *
   *  Counter |________|____|____|
   *            strm #   spwn ctr
   *
   *      | data[3]| data[2] |
   *
   *  Key |________|____|____|
   *       time stp| res|seed|
   *                 erv
   *                 ed
  */

  //! set the seed
  void set_seed(uint32_t const s){
    uint64_t tmp = data[2];
    tmp = tmp & ~SEED_MASK;
    data[2] = tmp | static_cast<uint64_t>(s);
    return;
  }

  //! set time step
  void set_time_step(uint64_t const s) {
    data[3] = s;
    return;
  }

  //! set stream number
  void set_stream_number(uint64_t const n){
    data[1] = n;
    return;
  }

  //! set spawn id
  void set_spawn_id(uint32_t const i){
    uint64_t tmp = data[0];
    // clear out any existing data, preserving the counter data
    tmp = tmp & ~SPAWN_MASK;
    data[0] = tmp | (static_cast<uint64_t>(i) << 32);
  }

  //! set the step counter
  void set_step_counter(uint32_t const c){
    uint64_t tmp = data[0];
    tmp = tmp & ~STEP_MASK;
    data[0] = tmp | static_cast<uint64_t>(c);
    return;
  }

  //! get the seed
  uint32_t get_seed() const { return get_seed_impl(data); }

  //! get step counter
  uint32_t get_step_counter() const { return get_step_counter_impl(data); }

  //! get time step
  uint64_t get_time_step() const { return get_time_step_impl(data); }

  //! get stream number
  uint64_t get_stream_number() const { return get_stream_number_impl(data); }

  //! get the spawn id
  uint32_t get_spawn_id() const { return get_spawn_id_impl(data); }

  bool can_spawn() const {
    return can_spawn_impl(data);
  }

  /**\brief Can this generator support spawning children? */
  // static bool can_spawn(ctr_type::value_type const *const rng_data) {
  //   // uint32_t const spawn_gen = get_spawn_gen_impl(rng_data);
  //   // uint32_t const max_gen = max_gens();
  //   // return spawn_gen < max_gen;
  //   return can_spawn_impl(rng_data);
  // }

  static constexpr uint32_t max_children_per_gen() { return 255; }

private:

  mutable ctr_type::value_type data[CBRNG_DATA_SIZE];

  mutable bool passed_max_ = false;

  //! Private copy constructor.
  Counter_RNG2(const Counter_RNG2 &);

  //! Private assignment operator.
  Counter_RNG2 &operator=(const Counter_RNG2 &);

  // private interface functions

  //! Initialize internal state from a seed and stream number.
  inline void initialize(const uint32_t seed, const uint64_t streamnum,
                         const uint64_t timestep);

  //! Spawn a new, independent generator from the provided state block.
  inline void _spawn(ctr_type::value_type const *const parent_data,
                     uint32_t const child_idx);

}; // class Counter_RNG2

//---------------------------------------------------------------------------//
// Implementation
//---------------------------------------------------------------------------//

//! Spawn a new, independent generator from this reference.
inline void Counter_RNG2_Ref::spawn(Counter_RNG2 &new_gen,
                                    uint32_t child_idx) const {
  new_gen._spawn(data.access(),child_idx);
}

//---------------------------------------------------------------------------//
//! Is this Counter_RNG2_Ref a reference to rng?
inline bool Counter_RNG2_Ref::is_alias_for(Counter_RNG2 const &rng) const {
  return rng.begin() == data.access();
}

//---------------------------------------------------------------------------//
//! \brief Initialize internal state from a seed, stream number, and time step.
//! \remark This should not be used in spawning!
inline void Counter_RNG2::initialize(const uint32_t seed,
                                     const uint64_t streamnum,
                                     const uint64_t timestep) {
  // If data[] is not an array of uint64_t, then this is likely broken!
  static_assert(std::is_same<ctr_type::value_type,uint64_t>::value,
    "Expected the ctr_type::value_type to be 64b unsigned int");

  set_step_counter(0);

  set_stream_number(streamnum);

  set_seed(seed);

  set_time_step(timestep);
}

//---------------------------------------------------------------------------//
/*! \brief Initialize this RNG's state by spawning from the provided parent
 * state block.
 *
 * \param parent_data: the data from the parent RNG.
 * \param idx: new child index, in the range [1,max_children_per_gen].
 *
 * \remark Check can_spawn() before spawning to ensure this generator is not
 * final generation.
 *
 * \remark Caller must provide an index of the child RNG.
 *
 * The child's spawn id, part of its RNG counter, is formed by ORing the
 * parent's id with a modifed version of the provided index. The provided
 * index is shifted left to move it one generation (8b) past the parent's gen.
 * The result is that the child will not conflict with any of the parent's
 * siblings.
 */
inline void Counter_RNG2::_spawn(ctr_type::value_type const *const parent_data,
                                 uint32_t idx) {
  Require(Counter_RNG2::can_spawn(parent_data));
  Require(idx > 0);
  Require(idx <= Counter_RNG2::max_children_per_gen());
  // Initialize this generator with the seed and stream number from the parent.
  uint32_t const seed = get_seed_impl(parent_data);
  uint64_t const streamnum = get_stream_number_impl(parent_data);
  uint64_t const time_step = get_time_step_impl(parent_data);

  this->initialize(seed, streamnum, time_step);

  uint32_t const parent_sid = get_spawn_id_impl(parent_data);
  uint32_t const parent_gen = get_spawn_gen_impl(parent_data);
  /* parent gen will not be greater than 3. Therefore, idx will be shifted
   * left at most 24b. */
  for(uint32_t i = 0; i < parent_gen; ++i){
    idx = idx << 8;
  }
  uint32_t const this_sid = parent_sid | idx;
  Check(this_sid > parent_sid);
  set_spawn_id(this_sid);
  return;
} // Counter_RNG2::_spawn(

} // end namespace rtt_rng

#endif

//---------------------------------------------------------------------------//
// end Counter_RNG2.hh
//---------------------------------------------------------------------------//
