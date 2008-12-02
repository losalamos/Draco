/*!
  \file    rng/LF_Gen.hh
  \author  Paul Henning
  \brief   Declaration of class LF_Gen
  \note    Copyright 2006 Los Alamos National Security, LLC.
  \version $Id$
*/

#ifndef LF_Gen_hh
#define LF_Gen_hh

#include <algorithm>
#include <ds++/Assert.hh>
#include <ds++/Data_Table.hh>
#include "rng/config.h"
#include "LFG.h"

namespace rtt_rng
{


class LF_Gen;

/*! This is a reference to an LF_Gen */
class LF_Gen_Ref
{
  public:
    LF_Gen_Ref(unsigned * const db, unsigned * const de)
        : data(db, de) 
    {
        Require(std::distance(db,de) == LFG_DATA_SIZE);
    }

    double ran() const
    {
        return lfg_gen_dbl(data.access()); 
    }

    //! Spawn a new, independent stream from this one.
    inline void spawn(LF_Gen& new_gen) const;

    //! Return the identifier for this stream
    unsigned get_num() const
    {
        return lfg_gennum(data.access());
    }


    inline bool is_alias_for(LF_Gen const &rng);


  private:
    mutable rtt_dsxx::Data_Table<unsigned> data;
};



/*! This holds the data for, and acts as the interface to, one random number
 * stream */
class LF_Gen
{
    friend class LF_Gen_Ref;
  public:
    typedef unsigned* iterator;
    typedef unsigned const * const_iterator;
  public:
    LF_Gen() 
    {
        Require(lfg_size() == LFG_DATA_SIZE);
    }


    LF_Gen(unsigned const seed, unsigned const streamnum)
    {
        // create a new Rnd object
        lfg_create_rng(streamnum, seed, begin(), end());
    }


    void finish_init() const
    {
        lfg_create_rng_part2(data, data + LFG_DATA_SIZE);
    }


    //! Return a random double
    double ran() const 
    { 
        return lfg_gen_dbl(data); 
    }


    //! Spawn a new, independent stream from this one.
    void spawn(LF_Gen& new_gen) const
    { 
        lfg_spawn_rng(data, new_gen.data, new_gen.data+LFG_DATA_SIZE); 
    }


    //! Return the identifier for this stream
    unsigned get_num() const
    {
        return lfg_gennum(data);
    }


    //! Return the size of the state
    unsigned size() const { return LFG_DATA_SIZE; }


    iterator begin() 
    { 
        return data; 
    }
    
    iterator end() 
    { 
        return data + LFG_DATA_SIZE; 
    }


    const_iterator begin() const 
    { 
        return data; 
    }

    const_iterator end() const 
    { 
        return data + LFG_DATA_SIZE; 
    }


    bool operator==(LF_Gen const & rhs) const 
    { 
        return std::equal(begin(), end(), rhs.begin()); 
    }

    LF_Gen_Ref ref() const
    {
        return LF_Gen_Ref(data, data+LFG_DATA_SIZE);
    }

    static unsigned size_bytes() { return LFG_DATA_SIZE*sizeof(unsigned); }
    
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



  private:
    mutable unsigned data[LFG_DATA_SIZE];
};

inline void LF_Gen_Ref::spawn(LF_Gen& new_gen) const
{ 
    lfg_spawn_rng(data.access(),  new_gen.data, new_gen.data+LFG_DATA_SIZE); 
}

inline bool LF_Gen_Ref::is_alias_for(LF_Gen const &rng)
{
    return rng.begin() == data.access();
}


}

#endif
