//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   linear/test/tstqrdcmp.cc
 * \author Kent Budge
 * \date   Mon Aug  9 13:39:20 2004
 * \brief  
 * \note   Copyright (C) 2016-2020 Triad National Security, LLC.
 *         All rights reserved.
 */
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

#include <iostream>
#include <vector>

#include "ds++/ScalarUnitTest.hh"
#include "ds++/Soft_Equivalence.hh"

#include "ds++/Release.hh"
#include "linear/qrdcmp.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_linear;

//----------------------------------------------------------------------------//
// TESTS
//----------------------------------------------------------------------------//

void tstqrdcmp(UnitTest &ut) {
  vector<double> A(4);
  A[0 + 2 * 0] = 2.;
  A[0 + 2 * 1] = 3.;
  A[1 + 2 * 0] = 1.;
  A[1 + 2 * 1] = 5.;

  vector<double> C, D;

  qrdcmp(A, 2, C, D);

  // Compute QR to verify
  double uj[2];
  uj[0] = A[0 + 2 * 0];
  uj[1] = A[1 + 2 * 0];
  double Qj[4];
  Qj[0 + 2 * 0] = 1 - uj[0] * uj[0] / C[0];
  Qj[0 + 2 * 1] = -uj[0] * uj[1] / C[0];
  Qj[1 + 2 * 0] = -uj[1] * uj[0] / C[0];
  Qj[1 + 2 * 1] = 1 - uj[1] * uj[1] / C[0];

  double QR[4];
  QR[0 + 2 * 0] = Qj[0 + 2 * 0] * D[0];
  QR[0 + 2 * 1] = Qj[0 + 2 * 0] * A[0 + 2 * 1] + Qj[0 + 2 * 1] * D[1];
  QR[1 + 2 * 0] = Qj[1 + 2 * 0] * D[0];
  QR[1 + 2 * 1] = Qj[1 + 2 * 0] * A[0 + 2 * 1] + Qj[1 + 2 * 1] * D[1];

  if (soft_equiv(QR[0 + 2 * 0], 2.0)) {
    ut.passes("0,0 is correct");
  } else {
    ut.failure("0,0 is NOT correct");
  }
  if (soft_equiv(QR[0 + 2 * 1], 3.0)) {
    ut.passes("0,1 is correct");
  } else {
    ut.failure("0,1 is NOT correct");
  }
  if (soft_equiv(QR[1 + 2 * 0], 1.0)) {
    ut.passes("1,0 is correct");
  } else {
    ut.failure("1,0 is NOT correct");
  }
  if (soft_equiv(QR[1 + 2 * 1], 5.0)) {
    ut.passes("1,1 is correct");
  } else {
    ut.failure("1,1 is NOT correct");
  }
}

//----------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  try {
    ScalarUnitTest ut(argc, argv, release);
    tstqrdcmp(ut);
  } catch (exception &err) {
    cout << "ERROR: While testing tstqrdcmp, " << err.what() << endl;
    return 1;
  } catch (...) {
    cout << "ERROR: While testing tstqrdcmp, an unknown exception was thrown."
         << endl;
    return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------//
// end of tstqrdcmp.cc
//----------------------------------------------------------------------------//
