//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   norms/test/tstL2norm.cc
 * \author Kent Budge
 * \date   Tue Sep 18 09:06:26 2007
 * \brief  Test the L2norm function template.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.  
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"
#include "norms/L2norm.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_norms;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstL2norm(UnitTest &ut) {
  unsigned const N = 2;
  unsigned const n = rtt_c4::nodes();

  vector<double> x(N, rtt_c4::node() + 1);

  double const norm = L2norm(x);

  if (soft_equiv(norm, sqrt(1. / 6 + n * (0.5 + n / 3.))))
    ut.passes("L2norm is correct");
  else
    ut.failure("L2norm is NOT correct");

  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstL2norm(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstL2norm.cc
//---------------------------------------------------------------------------//
