//----------------------------------*-C++-*----------------------------------//
// tstrn.cc
// Thomas M. Evans
// Wed Apr 29 15:11:32 1998
//---------------------------------------------------------------------------//
// @> Test of SPRNG random number class
//---------------------------------------------------------------------------//

#include "rng/Sprng.hh"
#include <iostream>

using namespace RNG;
using namespace std;

main()
{
    int *id1, *id2;
    int seed = 493875348;
    int num = 5;
    int **news, **ns;
    int spawn = 1;

    id1 = init_sprng(0, num, seed, 1);
    id2 = init_sprng(1, num, seed, 1);
    
    Sprng ran1(id1, 0);	
    Sprng ran2(ran1);

    for (int i = 0; i < 2; i++)
	cout << ran1.ran() << endl;

    for (int i = 2; i < 5; i++)
	cout << ran2.ran() << endl;

    ran1.print();
    ran2.print();

    cout << "**SPAWNING**" << endl;

    spawn = spawn_sprng(ran1.get_id(), 1, &news);

    for (int i = 0; i < 5; i++)
	cout << sprng(news[0]) << endl;

    Sprng ran3(news[0], 0);

    cout << endl;

    spawn = spawn_sprng(ran1.get_id(), 1, &ns);

    Sprng ran4(ns[0], 0);

    for (int i = 0; i < 5; i++)
	cout << ran4.ran() << endl;

    cout << news[0] << " " << ns[0] << endl;
}

//---------------------------------------------------------------------------//
//                              end of tstrn.cc
//---------------------------------------------------------------------------//
