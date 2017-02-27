//----------------------------------*-C++-*----------------------------------//
/*!
* \file   linear/test/tstfnorm.cc
* \author Kent Budge
* \date   Mon Aug  9 13:39:20 2004
* \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
*         All rights reserved.
*/
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "linear/fnorm.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void func(const vector<double> &x, vector<double> &fvec) {
  fvec.resize(2);
  fvec[0] = x[1] * sin(x[0]);
  fvec[1] = x[1] * cos(x[0]);
}

//---------------------------------------------------------------------------//
void tstfnorm(UnitTest &ut) {
  vector<double> x(2);
  vector<double> fvec;

  x[0] = 0.235;
  x[1] = 3.2;

  if (soft_equiv(fnorm(x, fvec, &func), 0.5 * 3.2 * 3.2)) {
    ut.passes("fnorm is correct");
  } else {
    ut.failure("fnorm is NOT correct");
  }
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstfnorm(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstfnorm.cc
//---------------------------------------------------------------------------//
