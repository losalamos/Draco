//----------------------------------*-C++-*----------------------------------//
// testrn2.cc
// Thomas M. Evans
// Thu Apr 30 13:33:20 1998
//---------------------------------------------------------------------------//
// @> Test of SPRNG and Sprng random number class
//---------------------------------------------------------------------------//

#include "rng/Sprng.hh"
#include "rng/Rnd_Control.hh"
#include <vector>
#include <iostream>

using namespace RNG;
using namespace std;

main()
{
    Rnd_Control rcon(493875348);
    vector<Sprng> bucket;

    {
	SP<Sprng> ran = rcon.get_rn();
	
	for (int i = 0; i < 2; i++)
	    cout << ran->ran() << endl;
    }

}

//---------------------------------------------------------------------------//
//                              end of testrn2.cc
//---------------------------------------------------------------------------//
