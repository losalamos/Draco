//----------------------------------*-C++-*----------------------------------//
// tstrn3.cc
// Thomas M. Evans
// Tue May  5 15:10:32 1998
//---------------------------------------------------------------------------//
// @> Simple random number test for cylinder sampling
//---------------------------------------------------------------------------//

#include "rng/Rnd_Control.hh"
#include "rng/Sprng.hh"
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;
using namespace RNG;

main()
{
    Rnd_Control rcon(493875348);
    Sprng ran = rcon.get_rn();

    double pi = 2.0 * asin(1.0);
    double area = 1;
    double r2 = sqrt(area/pi);
    double r22 = r2 * r2;
    double r1 = sqrt(r22/2);
    double r12 = 0.0;

    int binr1 = 0;
    int binr2 = 0;
    double r;

    int samples = 1000000;

    for (int i = 1; i <= samples; i++)
    {
	r = sqrt((r22 - r12) * ran.ran() + r12);
	r < r1 ? binr1++ : binr2++;
    }

    double r1t = static_cast<double>(binr1) / samples;
    double r2t = static_cast<double>(binr2) / samples;

    cout << "A1 Samples: " << r1t << endl;
    cout << "A2 Samples: " << r2t << endl;
}

//---------------------------------------------------------------------------//
//                              end of tstrn3.cc
//---------------------------------------------------------------------------//
