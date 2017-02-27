//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstto_string.cc
 * \author Kent Budge
 * \date   Fri Jul 25 08:49:48 2008
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/config.h"
#include "ds++/to_string.hh"
#include <cmath>   // M_PI on other platforms.
#include <cstdlib> // M_PI

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstto_string(UnitTest &ut) {
  string const pi = rtt_dsxx::to_string(M_PI, 20);

  if (soft_equiv(M_PI, atof(pi.c_str())))
    ut.passes("pi correctly written/read");
  else
    ut.failure("pi NOT correctly written/read");

  double const foo(2.11111111);
  unsigned int const p(23);
  // Must be careful to use rtt_dsxx::to_string and avoid std::to_string --
  // especially after 'using namespace std.'
  string s1(rtt_dsxx::to_string(foo, p));
  string s2(rtt_dsxx::to_string(foo));
  if (s1 == s2)
    ut.passes("double printed using default formatting.");
  else
    ut.failure("double printed with wrong format!");

  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    tstto_string(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstto_string.cc
//---------------------------------------------------------------------------//
