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

} // end namespace rtt_rng

//---------------------------------------------------------------------------//
//                              end of Sprng.cc
//---------------------------------------------------------------------------//
