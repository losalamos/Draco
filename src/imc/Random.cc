//----------------------------------*-C++-*----------------------------------//
// Random.cc
// Thomas M. Evans
// Thu Feb  5 13:34:55 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "imctest/Random.hh"
#include "imctest/Math.hh"
#include <iostream>
#include <fstream>

IMCSPACE

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// random number generator which returns a random deviate between 0.0 and 1.0
double Random::Ran()
{  
  // counter update
    count++;

  // variables and indices
    long mj, mk;
    int i,ii,k;

  // initialization
    if (idum < 0 || iff == 0)
    {
	iff = 1;
	mj = mseed - (idum < 0 ? -idum : idum);
	mj %= mbig;
	ma[55] = mj;
	mk = 1;
	for (i = 1; i <= 54; i++)
	{
	    ii = (21*i) % 55;
	    ma[ii] = mk;
	    mk = mj - mk;
	    if (mk < mz)
		mk += mbig;
	    mj = ma[ii];
	}
	for (k = 1; k <= 4; k++)
	{
	    for (i = 1; i <= 55; i++)
	    {
		ma[i] -= ma[1+(i+30) % 55];
		if (ma[i] < mz)
		    ma[i] += mbig;
	    }
	}
	inext = 0;
	inextp = 31;
	idum = 1;
    }

  // starting point after initialization
    if (++inext == 56) 
	inext = 1;
    if (++inextp == 56)
	inextp = 1;
    mj = ma[inext] - ma[inextp];
    if (mj < mz) mj += mbig;
    ma[inext] = mj;
    return static_cast<double>(mj) * fac;
}

//---------------------------------------------------------------------------//
// diagnostics and tests
//---------------------------------------------------------------------------//
// debug print
void Random::Print_values()
{
    using std::cout;
    using std::endl;
    using Global::operator<<;
    
    cout << "mbig:  " << mbig  << endl;
    cout << "mseed: " << mseed << endl;
    cout << "mz:    " << mz    << endl;
    cout << "fac:   " << fac   << endl;
    cout << "idum:  " << idum  << endl;
    cout << "count: " << count << endl;
    cout << "iff:   " << iff   << endl;
    cout << "ma:    " << endl;
    cout << ma;
}

//---------------------------------------------------------------------------//
// random number tests
double Random::Test_avg(int num)
{
    double avgcount = 0;
    for (int i = 1; i <= num; i++)
	avgcount += Ran();
    double average = avgcount / static_cast<double>(num);
    return average;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Random.cc
//---------------------------------------------------------------------------//
