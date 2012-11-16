//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Counter_RNG.hh
 * \author Peter Ahrens
 * \date   Fri Aug 3 16:53:23 2012
 * \brief  Declaration of class Counter_RNG.
 * \note   Copyright (C) 2012 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//

#ifndef Counter_RNG_hh
#define Counter_RNG_hh

#include "rng/config.h"

#ifdef _MSC_FULL_VER
// Engines have multiple copy constructors, quite legal C++, disable MSVC
// complaint.
#pragma warning (disable : 4521)
#endif

#if defined (__ICC)
// Suppress Intel's "unrecognized preprocessor directive" warning, triggered
// by use of #warning in Random123/features/sse.h.
#pragma warning disable 11
#endif

#define GNUC_VERSION (__GNUC__*10000 + __GNUC_MINOR__*100 + __GNUC_PATCHLEVEL__)
#if (GNUC_VERSION >= 40204) && !defined (__ICC) && !defined(NVCC)
// Suppress GCC's "unused parameter" warning, about lhs and rhs in sse.h, and
// an "unused local typedef" warning, from a pre-C++11 implementation of a
// static assertion in compilerfeatures.h.
#if (GNUC_VERSION >= 40600)
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

#include "Random123/threefry.h"

#if (GNUC_VERSION >= 40600)
// Restore GCC diagnostics to previous state.
#pragma GCC diagnostic pop
#endif

#include "Random123/u01.h"

#include <ds++/Data_Table.hh>
#include <algorithm>

namespace rtt_rng
{

// Forward declaration.
class Counter_RNG;

// Select a particular counter-based random number generator from Random123.
typedef r123::Threefry2x64 CBRNG;

#define CBRNG_DATA_SIZE 4

//===========================================================================//
/*!
 * \class Counter_RNG_Ref
 * \brief A reference to a Counter_RNG.
 *
 * Counter_RNG_Ref provides an interface to a counter-based random number
 * generator from the Random123 library from D. E. Shaw Research
 * (http://www.deshawresearch.com/resources_random123.html).  Unlike
 * Counter_RNG, Counter_RNG_Ref doesn't own its own RNG state (i.e., key and
 * counter); instead, it operates using a data block specified during
 * construction.
 */
//===========================================================================//
class Counter_RNG_Ref
{
  public:
    //! Constructor.  db and de specify the extents of an RNG state block.
    Counter_RNG_Ref(uint64_t* const db, uint64_t* const de)
        : data(db, de) 
    { Require(std::distance(db,de) * sizeof(uint64_t) ==
              sizeof(CBRNG::ctr_type) + sizeof(CBRNG::key_type)); }

    //! Return a random double in the open interval (0, 1)---i.e., excluding
    //! the endpoints.
    double ran() const
    {
        CBRNG rng;
        CBRNG::ctr_type ctr = {{data[0], data[1]}};
        CBRNG::key_type key = {{data[2], data[3]}};
        CBRNG::ctr_type result = rng(ctr, key);
        ctr.incr();
        std::copy(ctr.data(), ctr.data() + 2, data.access());
        return u01_open_open_64_53(result[0]);
    }

    //! Spawn a new, independent stream from this one.
    inline void spawn(Counter_RNG& new_gen) const;

    //! Return the identifier for this stream.
    uint64_t get_num() const { return data[2]; }

    //! Return a unique number for this stream and state.
    uint64_t get_unique_num() const { return data[2]; }

    //! Is this Counter_RNG_Ref a reference to rng?
    inline bool is_alias_for(Counter_RNG const &rng);

  private:
    mutable rtt_dsxx::Data_Table<uint64_t> data;
};


//===========================================================================//
/*!
 * \class Counter_RNG
 * \brief A counter-based random-number generator.
 *
 * Counter_RNG provides an interface to a counter-based random number
 * generator from the Random123 library from D. E. Shaw Research
 * (http://www.deshawresearch.com/resources_random123.html).
 */
//===========================================================================//
class Counter_RNG
{
  private:

    friend class Counter_RNG_Ref;
    mutable uint64_t data[4];

  public:

    typedef uint64_t* iterator;
    typedef uint64_t const * const_iterator;

    //! Default constructor.
    Counter_RNG() { Require(sizeof(data) ==
                            sizeof(CBRNG::ctr_type) +
                            sizeof(CBRNG::key_type)); }

    //! Constructor.
    Counter_RNG(uint64_t const seed, uint64_t const key_lo,
                uint64_t const key_hi)
    {
        data[0] = 0;
        data[1] = seed;
        data[2] = key_lo;
        data[3] = key_hi;
    }

    //! Create a new Counter_RNG from data.
    Counter_RNG(uint64_t* const _data)
    {
	std::copy(_data, _data + 4, data);
    }

    //! Return a random double in the interval (0, 1)---i.e., excluding the
    //! endpoints.
    double ran() const
    {
        CBRNG rng;
        CBRNG::ctr_type ctr = {{data[0], data[1]}};
        CBRNG::key_type key = {{data[2], data[3]}};
        CBRNG::ctr_type result = rng(ctr, key);
        ctr.incr();
        data[0] = ctr[0];
        data[1] = ctr[1];
        return u01_open_open_64_53(result[0]);
    }

    //! Spawn a new, independent stream from this one.
    void spawn(Counter_RNG& new_gen) const
    {
        std::copy(data, data + 4, new_gen.data);
        new_gen.data[0] = 0;   // Reset the lower counter in new stream.
        CBRNG::key_type key = {{new_gen.data[2], new_gen.data[3]}};
        key.incr();            // Increment the key in new stream.
        new_gen.data[2] = key[0];
        new_gen.data[3] = key[1];
    }

    //! Return the identifier for this stream.
    uint64_t get_num() const { return data[2]; }

    //! Return a unique number for this stream and state.
    uint64_t get_unique_num() const { return data[2]; }

    //! Return the size of the state.
    unsigned int size() const { return sizeof(data)/sizeof(uint64_t); }

    iterator begin() { return data; }
    
    iterator end() { return data + sizeof(data)/sizeof(uint64_t); }

    const_iterator begin() const { return data; }

    const_iterator end() const { return data + sizeof(data)/sizeof(uint64_t); }

    //! Test for equality.
    bool operator==(Counter_RNG const & rhs) const { 
        return std::equal(begin(), end(), rhs.begin()); }

    //! Return a Counter_RNG_Ref corresponding to this Counter_RNG.
    Counter_RNG_Ref ref() const {
        return Counter_RNG_Ref(data, data + sizeof(data)/sizeof(uint64_t)); }

    //! Return the size of this Counter_RNG in bytes.
    static unsigned int size_bytes() { return sizeof(data); }
    
  private:
    Counter_RNG(Counter_RNG const &);
};


//---------------------------------------------------------------------------//
// Implementation
//---------------------------------------------------------------------------//

// This implementation requires the full definition of Counter_RNG, so it must
// be placed here instead of in the Counter_RNG_Ref class.

inline void Counter_RNG_Ref::spawn(Counter_RNG& new_gen) const
{ 
    std::copy(data.begin(), data.begin() + 4, new_gen.data);
    new_gen.data[0] = 0;             // Reset the lower counter in new stream.
    CBRNG::key_type key = {{new_gen.data[2], new_gen.data[3]}};
    key.incr();                      // Increment the key in new stream.
    new_gen.data[2] = key[0];
    new_gen.data[3] = key[1];
}

inline bool Counter_RNG_Ref::is_alias_for(Counter_RNG const &rng)
{
    return rng.begin() == data.access();
}


} // end namespace rtt_rng

#endif
