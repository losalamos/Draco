//---------------------------------*-C++-*---------------------------------//
// random.hh
// Geoffrey Furnish
// 6 October 1992
//-------------------------------------------------------------------------//
// @> Decalre some functions for working with random numbers.  From
// @> Numerical Recipes.
//
// This file prototypes the functions in random.cc.  Eventually I'll
// need to go whole hog and make an actual random number class, but
// I'm in too much of a hurry to do that right now.  Besides, have
// to think about it some before implementing it.  Copied over from
// es1d. 
//-------------------------------------------------------------------------//

#ifndef __util_random_hh__
#define __util_random_hh__

float ran0(int *idum);
float ran1(int *idum);
float ran2(long *idum);
float ran3(int *idum);
float gasdev(int *idum);

#endif				// __util_random_hh__

//-------------------------------------------------------------------------//
//                              end of util/random.hh
//-------------------------------------------------------------------------//
