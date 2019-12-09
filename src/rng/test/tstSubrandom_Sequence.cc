//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/tstSubrandom_Sequence.cc
 * \author Kent Budge
 * \date   Thu Dec 22 14:16:45 2006
 * \brief  Test the Subrandom_Sequence class
 * \note   Copyright (C) 2006-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "rng/Halton_Sequence.hh"
#include "rng/Halton_Subrandom_Generator.hh"
#include "rng/LC_Subrandom_Generator.hh"
#include "rng/Sobol_Sequence.hh"
#include <fstream>
#include <iostream>

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_rng;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstSubrandom_Sequence(UnitTest &ut) {
  unsigned BASE = 1;
  unsigned COUNT = 25;

  Halton_Sequence seq(BASE);
  for (unsigned i = 0; i < COUNT; ++i) {
    cout << seq.shift() << endl;
  }

  double v1 = seq.lookahead();
  seq = Halton_Sequence(BASE, COUNT + 1);
  if (rtt_dsxx::soft_equiv(seq.lookahead(), v1)) {
    PASSMSG("correct restart");
  } else {
    FAILMSG("NOT correct restart");
  }

  {
    // Test Sobol_Sequence default constructor

    Sobol_Sequence sobol;
    // container should be empty
    std::vector<double> myValues(sobol.values());
    if (myValues.empty())
      FAILMSG(__LINE__);
    if (!soft_equiv(myValues[0], 0.5))
      FAILMSG(__LINE__);

    // Test Halton_Sequene default constructor

    Halton_Sequence halton;
    if (halton.base() != 0)
      FAILMSG(__LINE__);
    if (halton.count() != 0)
      FAILMSG(__LINE__);

    // Test LC_Subrandom_Generator

    LC_Subrandom_Generator lcsg;
    double value = lcsg.shift();
    if (!soft_equiv(value, 0.999741748906672))
      FAILMSG(__LINE__);
    lcsg.shift_vector();
  }

  Sobol_Sequence sobol(1);
  for (unsigned i = 0; i < COUNT; ++i) {
    cout << sobol.values()[0] << endl;
    sobol.shift();
  }

  // Sobol' sequence cannot presently be restarted.
  //     double v1 = seq.lookahead();
  //     seq = Sobol_Sequence(1, COUNT+1);
  //     if (seq.lookahead() == v1)
  //     {
  //         PASSMSG("correct restart");
  //     }
  //     else
  //     {
  //         FAILMSG("NOT correct restart");
  //     }

  // Integrate the unit sphere.
  double random_sum = 0.0;
  double subrandom_sum = 0.0;
  double sobol_sum = 0.0;
  Halton_Sequence xgen(0), ygen(1), zgen(2);
  Sobol_Sequence sgen(3);
  ofstream random_out("random.dat");
  ofstream subrandom_out("subrandom.dat");
  ofstream sobol_out("sobol.dat");
  Halton_Subrandom_Generator generator;
  unsigned const imax(100);
  for (unsigned i = 1; i < imax; ++i) {
    // random
    double x = 2 * (static_cast<double>(rand()) / RAND_MAX - 0.5);
    double y = 2 * (static_cast<double>(rand()) / RAND_MAX - 0.5);
    double z = 2 * (static_cast<double>(rand()) / RAND_MAX - 0.5);

    //         if (x*x+y*y+z*z<1.0) random_sum++;
    random_sum += (x + 1) * (x - 1) * (y + 1) * (y - 1) * (z + 1) * (z - 1);

    //         random_out << (8*random_sum/i - 4*PI/3) << '\n';
    random_out << (8 * random_sum / i + 64. / 27) << '\n';

    // subrandom
    x = 2 * (xgen.shift() - 0.5);
    y = 2 * (ygen.shift() - 0.5);
    z = 2 * (zgen.shift() - 0.5);

    //        if (x*x+y*y+z*z<1.0) subrandom_sum++;
    subrandom_sum += (x + 1) * (x - 1) * (y + 1) * (y - 1) * (z + 1) * (z - 1);

    // check generator
    if (!soft_equiv(x, 2 * (generator.shift() - 0.5)) ||
        !soft_equiv(y, 2 * (generator.shift() - 0.5)) ||
        !soft_equiv(z, 2 * (generator.shift() - 0.5))) {
      FAILMSG("Subrandom_Generator does not match Halton_Sequence");
      break;
    }
    generator.shift_vector();

    //         subrandom_out << (8*subrandom_sum/i - 4*PI/3) << '\n';
    subrandom_out << (8 * subrandom_sum / i + 64. / 27) << '\n';

    // Sobol
    x = 2 * (sgen.values()[0] - 0.5);
    y = 2 * (sgen.values()[1] - 0.5);
    z = 2 * (sgen.values()[2] - 0.5);
    sgen.shift();

    //        if (x*x+y*y+z*z<1.0) subrandom_sum++;
    sobol_sum += (x + 1) * (x - 1) * (y + 1) * (y - 1) * (z + 1) * (z - 1);

    //         subrandom_out << (8*subrandom_sum/i - 4*PI/3) << '\n';
    sobol_out << (8 * sobol_sum / i + 64. / 27) << '\n';
  }
  unsigned const myCount = generator.count();
  if (myCount == imax) {
    PASSMSG("Member function count returned the expected value (imax=100).");
  } else {
    FAILMSG(
        "Member function count did not return the expected value (imax=100).");
  }
  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstSubrandom_Sequence(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstSubrandom_Sequence.cc
//---------------------------------------------------------------------------//
