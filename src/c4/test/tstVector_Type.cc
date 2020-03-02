//----------------------------------*-C++-*-----------------------------------//
/*!
 * \file   c4/test/tstVector_Type.cc
 * \author Thomas M. Evans
 * \date   Tue Apr  2 15:57:11 2002
 * \brief  Ping Pong communication test.
 * \note   Copyright (C) 2002-2020 Triad National Security, LLC.
 *         All rights reserved. */
//----------------------------------------------------------------------------//

#include "c4/C4_Functions.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace std;
using namespace rtt_c4;

//----------------------------------------------------------------------------//
// TESTS
//----------------------------------------------------------------------------//

void test_simple(rtt_dsxx::UnitTest &ut) {
  // Create vector type for (1:2, 2:4) subset of a 4x6 double matrix.
  C4_Datatype data_type;
  int ierr = create_vector_type<double>(3, 2, 4, data_type);

  if (ierr == C4_SUCCESS)
    PASSMSG("created vector type successfully");
  else
    FAILMSG("did NOT create vector type successfully");

  // try sending a couple

  double array[4 * 6];

  unsigned proc = node();

  for (unsigned j = 0; j < 6; ++j) {
    for (unsigned i = 0; i < 4; ++i) {
      if (proc > 0)
        array[i + 4 * j] = 0;
      else
        array[i + 4 * j] = 10 * i + j;
    }
  }

  unsigned nproc = nodes();
  if (proc == 0) {
    for (unsigned p = 1; p < nproc; ++p)
      send_udt(array + (1 + 4 * 2), 1, p, data_type);
  } else {
    receive_udt(array + (1 + 4 * 2), 1, 0, data_type);
    for (unsigned j = 0; j < 6; ++j) {
      for (unsigned i = 0; i < 4; ++i) {
        if (i >= 1 && i <= 2 && j >= 2 && j <= 4) {
          if (!rtt_dsxx::soft_equiv(array[i + 4 * j], 10.0 * i + j)) {
            FAILMSG("did NOT transmit correct subarray");
            return;
          }
        } else {
          if (!rtt_dsxx::soft_equiv(array[i + 4 * j], 0.0)) {
            FAILMSG("did NOT transmit correct subarray");
            return;
          }
        }
      }
    }
    PASSMSG("transmitted correct subarray");
  }

  type_free(data_type);
}

//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_simple(ut);
  }
  UT_EPILOG(ut);
}

//----------------------------------------------------------------------------//
// end of tstVector_Type.cc
//----------------------------------------------------------------------------//
