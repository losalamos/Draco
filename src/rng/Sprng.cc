//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Sprng.cc
 * \author Thomas M. Evans
 * \date   Fri Jun 26 07:41:48 1998
 * \brief  \link rtt_rng::Sprng Sprng \endlink random number class 
 *         implementation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Sprng.hh"
#include <cmath>
#include <iostream>

namespace rtt_rng 
{

// stl components
using std::fabs;

int Sprng::packed_size = 0;

//---------------------------------------------------------------------------//
// diagnostic tests
//---------------------------------------------------------------------------//
// calculate the average for n deviates, should succeed

bool Sprng::avg_test(int n, double eps) const
{
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
    if (Sprng::packed_size > 0) return Sprng::packed_size;
    
    char *prng   = 0;
    int rng_size = pack_sprng(streamid->id, &prng);
    packed_size = rng_size + 2 * sizeof(int);
    Check (prng);

    return Sprng::packed_size;
    
}



} // end namespace rtt_rng

//---------------------------------------------------------------------------//
//                              end of Sprng.cc
//---------------------------------------------------------------------------//
