//----------------------------------*-C++-*----------------------------------//
/*!
 * \file    rng/LF_Gen.hh
 * \author  Paul Henning
 * \brief   Declaration of class LF_Gen
 * \note    Copyright (C) 2006-2013 Los Alamos National Security, LLC.
 *          All rights reserved.
 * \version $Id$
 */
//---------------------------------------------------------------------------//

#ifndef LF_Gen_hh
#define LF_Gen_hh

#include "LFG.h"
#include "ds++/Data_Table.hh"
#include "ds++/config.h"
#include <algorithm>

// A possible way to eliminate the warnings suppressed with this pragma is here: 
// http://www.windows-api.com/microsoft/VC-Language/30952961/a-solution-to-warning-c4251--class-needs-to-have-dllinterface.aspx
#if defined(MSVC)
#pragma warning (push)
// LF_Gen.hh(58): warning C4251: 'rtt_rng::LF_Gen_Ref::data' : class '
// rtt_dsxx::Data_Table<T>' needs to have dll-interface to be used by 
// clients of class 'rtt_rng::LF_Gen_Ref'
#pragma warning (disable:4251)
#endif

namespace rtt_rng
{

// forward declaration
class LF_Gen;

/*! This is a reference to an LF_Gen */
class DLL_PUBLIC LF_Gen_Ref
{
  public:
    LF_Gen_Ref(unsigned int* const db, unsigned int* const de)
        : data(db, de) 
    { Require(std::distance(db,de) == LFG_DATA_SIZE); }

    double ran() const { return lfg_gen_dbl(data.access()); }

    //! Spawn a new, independent stream from this one.
    inline void spawn(LF_Gen& new_gen) const;

    //! Return the identifier for this stream
    unsigned int get_num() const { return lfg_gennum(data.access()); }

    //! Return a unique number for this stream and state
    unsigned int get_unique_num() const { return lfg_unique_num(data.access()); }

    inline bool is_alias_for(LF_Gen const &rng);

  private:

    mutable rtt_dsxx::Data_Table<unsigned int> data;
};

//===========================================================================//
/*!
 * \class LF_Gen
 * \brief This holds the data for, and acts as the interface to, one random
 *        number stream
 */
//===========================================================================//
class DLL_PUBLIC LF_Gen
{
  private:

    friend class LF_Gen_Ref;
    mutable unsigned int data[LFG_DATA_SIZE];

  public:

    typedef unsigned int* iterator;
    typedef unsigned int const * const_iterator;

    //! Default constructor.
    LF_Gen() { Require(lfg_size() == LFG_DATA_SIZE); }

    LF_Gen(unsigned int const seed, unsigned int const streamnum)
    {
        // create a new Rnd object
        lfg_create_rng(streamnum, seed, begin(), end());
    }

    LF_Gen(unsigned int* const _data)
    {
#if LFG_PARAM_SET == 1
        Require( _data[LFG_DATA_SIZE-4] < 17 );
#elif LFG_PARAM_SET == 2
        Require( _data[LFG_DATA_SIZE-4] < 31 );
#endif
	// create a new Rnd object from data
	std::copy (_data, _data + LFG_DATA_SIZE, data);
    }

    void finish_init() const {
        lfg_create_rng_part2(data /*, data + LFG_DATA_SIZE */); }

    //! Return a random double
    double ran() const { return lfg_gen_dbl(data); }

    //! Spawn a new, independent stream from this one.
    void spawn(LF_Gen& new_gen) const { 
        lfg_spawn_rng(data, new_gen.data, new_gen.data+LFG_DATA_SIZE); }

    //! Return the identifier for this stream
    unsigned int get_num() const { return lfg_gennum(data); }

    //! Return a unique number for this stream and state
    unsigned int get_unique_num() const { return lfg_unique_num(data); }

    //! Return the size of the state
    unsigned int size() const { return LFG_DATA_SIZE; }

    iterator begin() { return data; }
    
    iterator end() { return data + LFG_DATA_SIZE; }

    const_iterator begin() const { return data; }

    const_iterator end() const { return data + LFG_DATA_SIZE; }

    bool operator==(LF_Gen const & rhs) const { 
        return std::equal(begin(), end(), rhs.begin()); }

    LF_Gen_Ref ref() const {
        return LF_Gen_Ref(data, data+LFG_DATA_SIZE); }

    static unsigned int size_bytes() {
        return LFG_DATA_SIZE*sizeof(unsigned int); }
    
#if 0
    // Copying RNG streams shouldn't be done lightly!
    inline LF_Gen& operator=(LF_Gen const &src)
    {
        if(&src != this)
            std::memcpy(data, src.data, size_bytes());
        return *this;
    }
#endif

  private:
    LF_Gen(LF_Gen const &);

};

//---------------------------------------------------------------------------//
// Implementation
//---------------------------------------------------------------------------//

// This implementation requires the full definition of LF_Gen, so it must be
// placed here instead of in the LF_Gen_Ref class.

inline void LF_Gen_Ref::spawn(LF_Gen& new_gen) const
{ 
    lfg_spawn_rng(data.access(), new_gen.data, new_gen.data+LFG_DATA_SIZE); 
}

inline bool LF_Gen_Ref::is_alias_for(LF_Gen const &rng)
{
    return rng.begin() == data.access();
}

} // end namespace rtt_rng

#if defined(MSVC)
#   pragma warning (pop)
#endif


#endif

//---------------------------------------------------------------------------//
// end of LF_Gen.hh
//---------------------------------------------------------------------------//