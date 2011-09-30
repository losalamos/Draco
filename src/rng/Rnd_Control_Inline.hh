//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Rnd_Control_Inline.hh
 * \author Paul Henning
 * \brief  Rnd_Control header file.
 * \note   Copyright (C) 2008-2011 LANS, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_rng_Rnd_Control_Inline_hh
#define rtt_rng_Rnd_Control_Inline_hh

//#include "ds++/Assert.hh"
//#include "rng/config.h"

// header file for SPRNG package
//#include "LFG.h"
#include "LF_Gen.hh"

namespace rtt_rng 
{

/*! All this does is to hold a seed and a stream number.  Everything else
 *  that the traditional Rnd_Control did was moved to LF_Gen */

class Rnd_Control 
{
  private:
    // seed for initialization of random number streams
    unsigned d_seed;

    // total number of streams allowed
    const unsigned d_number;

    // number of current stream
    unsigned d_streamnum;

  public:
    // Constructor.
    Rnd_Control(int seed, int max_streams = 1000000000, int sn = 0, int p = 1);


    //! Query for the current random number stream index.
    int get_num() const { return d_streamnum; }

    //! Set (reset) the random number stream index.
    void set_num(const int num) { d_streamnum = num; }

    //! Query size of a packed random number state.
    int get_size() const { return (LFG_DATA_SIZE+2)*sizeof(unsigned); }

    //! Get the seed value used to initialize the SPRNG library.
    int get_seed() const { return d_seed; }

    //! Return the total number of current streams set.
    int get_number() const { return d_number; }

    inline void initialize(const unsigned snum, LF_Gen&);
    inline void initialize(LF_Gen&);
    inline void half_initialize(LF_Gen&);
};


// ---------------------------------------------------------------------------


//! Update the stream number and initialize the LF_Gen
inline void 
Rnd_Control::initialize(const unsigned snum, LF_Gen& lf)
{
    Require (snum <= d_number);

    // reset streamnum
    d_streamnum = snum;

    // create a new Rnd object
    lfg_create_rng(d_streamnum, d_seed, lf.begin(), lf.end());

    // advance the counter
    d_streamnum++;
}


//! Initialize the LF_Gen with the next stream number
inline void 
Rnd_Control::initialize(LF_Gen& lf)
{
    Require(d_streamnum <= d_number);

    // create a new Rnd object
    lfg_create_rng(d_streamnum, d_seed, lf.begin(), lf.end());

    // advance the counter
    d_streamnum++;

}


//! Do half an initialization... it is your responsibility to finish!
inline void 
Rnd_Control::half_initialize(LF_Gen& lf)
{
    Require(d_streamnum <= d_number);

    // create a (partially initialized) new Rnd object
    lfg_create_rng_part1(d_streamnum, d_seed, lf.begin(), lf.end());

    // advance the counter
    d_streamnum++;
}



} // end namespace rtt_rng

#endif                          // rtt_rng_Rnd_Control_hh

//---------------------------------------------------------------------------//
//                              end of rng/Rnd_Control.hh
//---------------------------------------------------------------------------//
