//----------------------------------*-C++-*----------------------------------//
// tstSprng.cc
// Thomas M. Evans
// Wed Apr 29 15:11:32 1998
// $Id$
//---------------------------------------------------------------------------//
// @> Test of SPRNG random number class
//---------------------------------------------------------------------------//

#include "RNG_Test.hh"
#include "../Random.hh"
#include "../Release.hh"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;
using rtt_rng::Sprng;

// global stuff
int seed = 493875348;

// passing condition
bool passed = true;
#define ITFAILS passed = rtt_rng_test::fail(__LINE__);

//---------------------------------------------------------------------------//

void ran_test()
{
    // data for random numbers
    int *id1, *id2, *id3;
    int num = 5;

    // get two random number seeds
    id1 = init_sprng(0, num, seed, 1);
    id2 = init_sprng(1, num, seed, 1);
    id3 = init_sprng(0, num, seed, 1);
    
    Sprng ran1(id1, 0);	
    Sprng ran2(id2, 1);

    if (!ran1.avg_test(10000, .01)) ITFAILS;
    if (!ran2.avg_test(10000, .01)) ITFAILS;

    Sprng ran3(id3, 0);

    if (!ran3.avg_test(10000, .01)) ITFAILS;

    // now check that 1, 3 give equal streams
    double eps = .00001;
    for (int i = 0; i < 10; i++)
	if (fabs(ran1.ran() - ran3.ran()) > eps) ITFAILS;
}

//---------------------------------------------------------------------------//

void sprng_test()
{
    int num  = 5;

    // make two sprng states
    int *idr = init_sprng(0, num, seed, 1);
    int *id1 = init_sprng(0, num, seed, 1);

    // now make some sprngs
    Sprng ranr(idr, 0);
    Sprng ran1(id1, 0);
    Sprng ran2(ran1);
    Sprng ran3(ran1);

    // get some reference numbers
    vector<double> ref(80);
    for (int i = 0; i < 80; i++)
	ref[i] = ranr.ran();

    // now check these sprngs that ALL access the same random number state
    double eps = .00001;
    for (int i = 0; i < 20; i++)
	if (fabs(ran1.ran() - ref[i]) > eps) ITFAILS;
    for (int i = 20; i < 40; i++)
	if (fabs(ran2.ran() - ref[i]) > eps) ITFAILS;
    for (int i = 40; i < 60; i++)
	if (fabs(ran3.ran() - ref[i]) > eps) ITFAILS;

    // now check the ids
    if (ran1.get_id() != id1) ITFAILS;
    if (ran2.get_id() != id1) ITFAILS;
    if (ran3.get_id() != id1) ITFAILS;

    // now assignment
    ranr = ran2;
    for (int i = 60; i < 80; i++)
	if (fabs(ranr.ran() - ref[i]) > eps) ITFAILS;
    if (ranr.get_id() != id1) ITFAILS;
}

//---------------------------------------------------------------------------//

void spawn_test()
{
//     int **news;
//     int spawn = 1;

//     spawn = spawn_sprng(ran1.get_id(), 1, &news);
    
//     for (int i = 0; i < 5; i++)
// 	cout << sprng(news[0]) << endl;
    
//     free_sprng(news[0]);
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    // version tag
    for (int arg = 1; arg < argc; arg++)
	if (string(argv[arg]) == "--version")
	{
	    cout << argv[0] << ": version " << rtt_rng::release() << endl; 
	    return 0;
	}

    // tests
    ran_test();
    sprng_test();

    // status of test
    cout << endl;
    cout <<     "*********************************" << endl;
    if (passed) 
    {
        cout << "**** Sprng Self Test: PASSED ****" << endl;
    }
    cout <<     "*********************************" << endl;
    cout << endl;
    
    cout << "Done testing Sprng." << endl;
}
//---------------------------------------------------------------------------//
//                              end of tstSprng.cc
//---------------------------------------------------------------------------//
