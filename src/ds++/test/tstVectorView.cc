//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstVectorView.cc
 * \author Kelly Thompson <kgt@lanl.gov>
 * \date   Sunday, Dec 15, 2019, 17:50 pm
 * \note   Copyright (C) 2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/VectorView.hh"
#include <limits>
#include <sstream>

using namespace std;
using namespace rtt_dsxx;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
void test_int_data(ScalarUnitTest &ut) {

  // 120 length vector
  vector<int> myvec(4 * 5 * 6, 0);
  for (size_t i = 1; i < myvec.size(); ++i) // make value == index
    myvec[i] = myvec[i - 1] + 1;

  // test entries for 3-D view.
  VectorView<int, 4, 5, 6> myvecview(
      myvec); // recast 1D as 3D with these lengths.
  FAIL_IF_NOT(myvecview(2, 0, 0) == 2);
  FAIL_IF_NOT(myvecview(0, 2, 0) == 8);
  FAIL_IF_NOT(myvecview(0, 0, 4) == 80);
  FAIL_IF_NOT(myvecview(2, 3, 4) == 94);
  FAIL_IF_NOT(myvecview(3, 4, 5) == 119);
}

//----------------------------------------------------------------------------//
void test_double_data(ScalarUnitTest &ut) {

  // 120 length vector
  vector<double> myvec(4 * 5 * 6, 0.0);
  for (size_t i = 1; i < myvec.size(); ++i) // make value == index
    myvec[i] = myvec[i - 1] + 1.0001;

  // test entries for 3-D view.
  VectorView<double, 4, 5, 6> myvecview(
      myvec); // recast 1D as 3D with these lengths.
  FAIL_IF_NOT(soft_equiv(myvecview(2, 0, 0), 2.0002));
  FAIL_IF_NOT(soft_equiv(myvecview(0, 2, 0), 8.0008));
  FAIL_IF_NOT(soft_equiv(myvecview(0, 0, 4), 80.0080));
  FAIL_IF_NOT(soft_equiv(myvecview(2, 3, 4), 94.0094));
  FAIL_IF_NOT(soft_equiv(myvecview(3, 4, 5), 119.0119));
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, release);
  try {
    test_int_data(ut);
    test_double_data(ut);
    ut.passes("Just Because.");
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstVectorView.cc
//---------------------------------------------------------------------------//
