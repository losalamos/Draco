//----------------------------------*-C++-*----------------------------------//
// tstrn4.cc
// Thomas M. Evans
// Wed Jun  3 13:38:02 1998
//---------------------------------------------------------------------------//
// @> Simple random number test for tilting
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
  // cell defs
    double m  = 2;
    double b  = 5;
    double x1 = 10;
    double x2 = 20;
    double a  = m/2. * (x2*x2 - x1*x1) + b * (x2 - x1);
    double p1 = 25. * (x2 - x1) / a;
    double c  = 1 / (m/2 * (x2*x2 - x1*x1) + b * (x2 - x1));

  // make random number object
    Rnd_Control rcon(493875348);
    Sprng ran = rcon.get_rn();

  // set bins
    vector<int> bin1(11), bin2(11);

  // sample
    for (int i = 1; i <= 1000000; i++)
    {
	double s1, s2;

      // sample the first way
	if (ran.ran() < p1)
	    s1 = x1 + (x2 - x1) * ran.ran();
	else
	    s1 = x1 + (x2 - x1) * sqrt(ran.ran());

      // sample the second way
	s2 = (-c*b + sqrt((c*b)*(c*b) + 4*(c*m/2)*(ran.ran() + c*m*x1*x1/2 + 
						  c*b*x1))) / (c*m);

      // bin them up
	if (s1 >= 10 && s1 < 11) bin1[1]++;
	if (s1 >= 11 && s1 < 12) bin1[2]++;
	if (s1 >= 12 && s1 < 13) bin1[3]++;
	if (s1 >= 13 && s1 < 14) bin1[4]++;	
	if (s1 >= 14 && s1 < 15) bin1[5]++;
	if (s1 >= 15 && s1 < 16) bin1[6]++;
	if (s1 >= 16 && s1 < 17) bin1[7]++;
	if (s1 >= 17 && s1 < 18) bin1[8]++;	
	if (s1 >= 18 && s1 < 19) bin1[9]++;
	if (s1 >= 19 && s1 < 20) bin1[10]++;

	if (s2 >= 10 && s2 < 11) bin2[1]++;
	if (s2 >= 11 && s2 < 12) bin2[2]++;
	if (s2 >= 12 && s2 < 13) bin2[3]++;
	if (s2 >= 13 && s2 < 14) bin2[4]++;	
	if (s2 >= 14 && s2 < 15) bin2[5]++;
	if (s2 >= 15 && s2 < 16) bin2[6]++;
	if (s2 >= 16 && s2 < 17) bin2[7]++;
	if (s2 >= 17 && s2 < 18) bin2[8]++;	
	if (s2 >= 18 && s2 < 19) bin2[9]++;
	if (s2 >= 19 && s2 < 20) bin2[10]++;
    }

    for (int i = 1; i <= 10; i++)
    {
	cout << "Method 1: " << i << " " << static_cast<double>(bin1[i]) /
	    1000000.0 << endl;
	cout << "Method 2: " << i << " " << static_cast<double>(bin2[i]) /
	    1000000.0 << endl;
	cout << endl;
    }
}
//---------------------------------------------------------------------------//
//                              end of tstrn4.cc
//---------------------------------------------------------------------------//
