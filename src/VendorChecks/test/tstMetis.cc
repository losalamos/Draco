//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   VendorChecks/test/tstMetis.cc
 * \date   Wednesday, May 11, 2016, 12:01 pm
 * \brief  Attempt to link to libmetis and run a simple problem.
 * \note   Copyright (C) 2016, Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <metis.h>
#include <vector>

// Original provided by Erik Zenker
// https://gist.github.com/erikzenker/c4dc42c8d5a8c1cd3e5a

void test_metis(rtt_dsxx::UnitTest &ut) {
  int nVertices = 6;
  int nWeights = 1;
  int nParts = 2;

  int objval(0);
  std::vector<int> part(nVertices, 0);

  // Indexes of starting points in adjacent array
  int xadj[] = {0, 2, 5, 7, 9, 12, 14};

  // Adjacent vertices in consecutive index order
  int adjncy[] = {1, 3, 0, 4, 2, 1, 5, 0, 4, 3, 1, 5, 4, 2};

  // Weights of vertices
  // if all weights are equal then can be set to NULL
  std::vector<int> vwgt(nVertices * nWeights, 0);

  // Partition a graph into k parts using either multilevel recursive
  // bisection or multilevel k-way partitioning.
  int ret =
      METIS_PartGraphKway(&nVertices, &nWeights, xadj, adjncy, NULL, NULL, NULL,
                          &nParts, NULL, NULL, NULL, &objval, &part[0]);

  if (ret == METIS_OK)
    PASSMSG("Successfully called METIS_PartGraphKway().");
  else
    FAILMSG("Call to METIS_PartGraphKway() failed.");

  int expectedResult[] = {1, 1, 0, 1, 0, 0};
  std::vector<int> vExpectedResult(expectedResult, expectedResult + 6);
  if (part == vExpectedResult)
    PASSMSG("Metis returned the expected result.");
  else
    FAILMSG("Metis failed to return the expected result.");

  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_metis(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstMetis.cc
//---------------------------------------------------------------------------//
