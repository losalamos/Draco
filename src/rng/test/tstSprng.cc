//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/tstSprng.cc
 * \author Thomas M. Evans
 * \date   Mon Dec 17 16:05:26 2001
 * \brief  Sprng testing.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "rng_test.hh"
#include "../Release.hh"
#include "../Random.hh"
#include "ds++/Assert.hh"
#include "ds++/Soft_Equivalence.hh"

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
using rtt_rng::Sprng;

// global stuff
int seed = 493875348;

using namespace std;

//---------------------------------------------------------------------------//
// TESTS
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

    if (rtt_rng_test::passed)
	PASSMSG("Simple random number access test passes.");
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

    if (rtt_rng_test::passed)
	PASSMSG("Simple Sprng object test passes.");
}

//---------------------------------------------------------------------------//

void pack_test()
{
    int num = 5;

    vector<double> ref(80);
    vector<char>   packed;

    {
	// make a sprng state
	int *id1 = init_sprng(0, num, seed, 1);
	int *idr = init_sprng(0, num, seed, 1);
	
	// now make some sprngs
	Sprng ran1(id1, 0);
	Sprng ranr(idr, 0);

	// get some reference numbers
	for (int i = 0; i < 80; i++)
	    ref[i] = ranr.ran();

	// get 40 numbers from the non-ref
	for (int i = 0; i < 40; i++)
	    ran1.ran();

	// pack up the sprng
	packed = ran1.pack();
    }

    // now check it
    {
	Sprng uran(packed);

	// now check the first 40 Unpacked ran numbers
	double r   = 0;
	double rf  = 0;
	for (int i = 0; i < 40; i++)
	{
	    r  = uran.ran();
	    rf = ref[i+40];
 
	    if (!rtt_dsxx::soft_equiv(r,rf)) ITFAILS;
	}
    }

    if (rtt_rng_test::passed)
	PASSMSG("Packing/Unpacking ok.");
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
	    cout << argv[0] << ": version " << rtt_rng::release() 
		 << endl;
	    return 0;
	}

    try
    {
	// >>> UNIT TESTS
	ran_test();
	sprng_test();
	pack_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
	cout << "While testing tstSprng, " << ass.what()
	     << endl;
	return 1;
    }

    // status of test
    cout << endl;
    cout <<     "*********************************************" << endl;
    if (rtt_rng_test::passed) 
    {
        cout << "**** tstSprng Test: PASSED" 
	     << endl;
    }
    cout <<     "*********************************************" << endl;
    cout << endl;

    cout << "Done testing tstSprng." << endl;
}   

//---------------------------------------------------------------------------//
//                        end of tstSprng.cc
//---------------------------------------------------------------------------//
