//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstRange_finder.cc
 * \author Mike Buksas
 * \date   Thu Feb  6 12:43:22 2003
 * \note   Copyright (C) 2016-2018 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Range_Finder.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void test_range_finder_left(UnitTest &ut) {

  double v[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

  vector<double> values(v, v + 10);

  int index = Range_finder_left(v, v + 10, 1.5);
  if (index != 1)
    ut.failure("test FAILS");

  index = Range_finder_left(values.begin(), values.end(), 2.5);
  if (index != 2)
    ut.failure("test FAILS");

  // Check for equality at all values:
  for (int i = 0; i < 10; ++i) {
    index = Range_finder_left(v, v + 10, static_cast<double>(i));
    if (index != i)
      ut.failure("test FAILS");
  }

  // For equality with the last value, we should get n-1 with end catching:
  index = Range_finder_left_catch_end(v, v + 10, 9.0);
  if (index != 8)
    ut.failure("test FAILS");

  index = Range_finder_catch_end(v, v + 10, 9.0, LEFT);
  if (index != 8)
    ut.failure("test FAILS");

  //     index = Range_finder_left(v,v+10, 42.69);
  //     if (index != -1) ut.failure("test FAILS");

  //     index = Range_finder_left(v+5,v+10, 1.0);
  //     if (index != -1) ut.failure("test FAILS");

  double rv[10] = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0};

  vector<double> rvalues(rv, rv + 10);

  index = Range_finder(rvalues.rbegin(), rvalues.rend(), 5.5, LEFT);
  if (index != 5)
    ut.failure("test FAILS");

  index = Range_finder_left(rvalues.rbegin(), rvalues.rend(), 5.0);
  if (index != 5)
    ut.failure("test FAILS");

  //     index = Range_finder_left(rvalues.rbegin(), rvalues.rend(), 10.12);
  //     if (index != -1) ut.failure("test FAILS");

  if (validate(std::pair<double *, double *>(rv + 0, rv + 0), rv + 0,
               rv + 10)) {
    ut.failure("validate FAILED to catch out of range result");
  } else {
    ut.passes("validate caught out of range result");
  }
  if (validate(std::pair<double *, double *>(rv + 10, rv + 10), rv + 0,
               rv + 10)) {
    ut.failure("validate FAILED to catch out of range result");
  } else {
    ut.passes("validate caught out of range result");
  }
  return;
}

void test_range_finder_right(UnitTest &ut) {

  double v[10] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

  int index;

  // Check for equality at all values. Note that 0 comes back as interval
  // -1 (e.g. out of range).
  for (int i = 1; i < 10; ++i) {
    index = Range_finder(v, v + 10, static_cast<double>(i), RIGHT);
    if (index != i - 1)
      ut.failure("test FAILS");
  }

  index = Range_finder_right_catch_end(v, v + 10, 0.0);
  if (index != 0)
    ut.failure("test FAILS");

  index = Range_finder_catch_end(v, v + 10, 0.0, RIGHT);
  if (index != 0)
    ut.failure("test FAILS");
  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    test_range_finder_left(ut);
    test_range_finder_right(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstRange_finder.cc
//---------------------------------------------------------------------------//
