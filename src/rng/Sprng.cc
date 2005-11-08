//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Sprng.cc
 * \author Thomas M. Evans
 * \date   Fri Jun 26 07:41:48 1998
 * \brief  Sprng random number class implementation file.
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <cmath>
#include <iostream>
#include <cstdlib>
#include "ds++/Packing_Utils.hh"
#include "Sprng.hh"

namespace rtt_rng 
{

// stl components
using std::fabs;

int Sprng::packed_size = 0;

//---------------------------------------------------------------------------//
/*!
 * \brief Unpacking constructor.
 */
Sprng::Sprng(const std::vector<char> &packed)
    : streamid(0), streamnum(0)
{
    Require (packed.size() >= 2 * sizeof(int));

    // make an unpacker
    rtt_dsxx::Unpacker u;
    
    // set the buffer
    u.set_buffer(packed.size(), &packed[0]);

    // unpack the stream num and size of state
    int rng_size = 0;
    u >> streamnum >> rng_size;
    Check (streamnum >= 0);
    Check (rng_size >= 0);

    // unpack the random stream state
    char *prng = new char[rng_size];
    for (int i = 0; i < rng_size; i++)
	u >> prng[i];

    // now rebuild the sprng object
    int *rnid = unpack_sprng(prng);

    // now make a new streamid
    streamid = new SprngValue(rnid);

    // reclaim memory
    delete [] prng;

    Ensure (u.get_ptr() == &packed[0] + packed.size());
    Ensure (streamid);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Pack a Sprng object into a vector<char>.
 */
std::vector<char> Sprng::pack() const
{
    Require (streamid);
    
    // make a packer
    rtt_dsxx::Packer p;

    // first pack the random number object and determine the size
    char *prng   = 0;
    int rng_size = pack_sprng(streamid->id, &prng);
    int size     = rng_size + 2 * sizeof(int);
    Check (prng);

    // now set the buffer
    std::vector<char> packed(size);
    p.set_buffer(size, &packed[0]);

    // pack the stream number and rng size
    p << streamnum << rng_size;

    // pack the stream state
    for (int i = 0; i < rng_size; i++)
	p << prng[i];

    // free the prng buffer
    std::free(prng);

    Ensure (p.get_ptr() == &packed[0] + size);

    return packed;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment operator.
 */
Sprng& Sprng::operator=(const Sprng &rhs)
{
    Require(rhs.streamid);

    // check to see if the values are the same
    if (streamid == rhs.streamid && streamnum == rhs.streamnum)
	return *this;

    // destroy this' value if it was the last random number in a particular
    // stream 
    if (streamid && --streamid->refcount == 0)
	delete streamid;

    // do assignment
    streamid  = rhs.streamid;
    streamnum = rhs.streamnum; 
    ++streamid->refcount;
    return *this;
}

//---------------------------------------------------------------------------//
// diagnostic tests
//---------------------------------------------------------------------------//
// calculate the average for n deviates, should succeed

bool Sprng::avg_test(int n, double eps) const
{
    Require(streamid);
    double avg = 0.0;
    for (int i = 1; i <= n; i++)
	avg += ran();
    double result = avg / n - 0.5;
    if (fabs(result) >= eps)
	return false;
    return true;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the size of the packed character stream.
 * 
 * \return The size
 */
int Sprng::get_size() const
{
    Require(streamid);
    if (Sprng::packed_size > 0) return Sprng::packed_size;
    
    char *prng   = 0;
    int rng_size = pack_sprng(streamid->id, &prng);
    packed_size  = rng_size + 2 * sizeof(int);

    // clear memory
    std::free(prng);

    return Sprng::packed_size;
}

} // end namespace rtt_rng

//---------------------------------------------------------------------------//
//                              end of Sprng.cc
//---------------------------------------------------------------------------//
