//----------------------------------*-C++-*----------------------------------//
// Random.hh
// Thomas M. Evans
// Thu Feb  5 13:34:55 1998
//---------------------------------------------------------------------------//
// @> Random number class header
//---------------------------------------------------------------------------//

#ifndef __imctest_Random_hh__
#define __imctest_Random_hh__

//===========================================================================//
// class Random - 
//
// Purpose : Random number class for use with IMC and MC applications;
//           origin, RAN3 function, pg.283 Numerical Recipes in C
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imctest/Names.hh"
#include <vector>
#include <algorithm>
#include <iostream>

IMCSPACE

using std::vector;
using std::fill;
using std::ostream;

class Random
{
private:
  // generator const variables
    const long mbig;
    const long mseed;
    const int mz;
    const double fac;
  // variables used by Ran()
    int inext, inextp;
    vector<long> ma;
    int iff;
  // seed value, negative integer initializes the sequence
    long idum;
  // counter
    long count;
  // original seed value
    const long seed;

public:
  // must give seed to Random number object
    explicit inline Random(long);
 
  // overloaded assignment operators
    inline const Random& operator=(const Random &rhs);

  // get Random number function
    double ran();

  // determine the number of random numbers used in this object
    long get_count() { return count; }

  // random number diagnostics
    long get_seed() { return seed; }
    void print_values(); 
    double test_avg(int num);
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

// overloaded stream operator
inline ostream& operator<<(ostream &output, Random &object)
{
    output << object.ran();
    return output;
}

//---------------------------------------------------------------------------//

inline const Random& Random::operator=(const Random &rhs)
{
  // overloaded assignment operator
    inext  = rhs.inext;
    inextp = rhs.inextp;
    ma     = rhs.ma;
    iff    = rhs.iff;
    idum   = rhs.idum;
    count  = rhs.count;

    return *this;
}

//---------------------------------------------------------------------------//
// inline functions for Random
//---------------------------------------------------------------------------//

// constructor
inline Random::Random(long idum_)
    : mbig(1000000000), mseed(161803398), mz(0), fac(1.0/mbig),
      inext(0), inextp(0), ma(56), iff(0), idum(idum_), count(0),
      seed(idum_)
{ 
    fill(ma.begin(), ma.end(), 0); 
}

CSPACE

#endif                          // __imctest_Random_hh__

//---------------------------------------------------------------------------//
//                              end of imctest/Random.hh
//---------------------------------------------------------------------------//
