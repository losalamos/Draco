//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   VendorChecks/test/tstMetis.cc
 * \date   Wednesday, May 11, 2016, 12:01 pm
 * \brief  Attempt to link to libmetis and run a simple problem.
 * \note   Copyright (C) 2016-2019, Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include <metis.h>
#include <vector>

// Original provided by Erik Zenker
// https://gist.github.com/erikzenker/c4dc42c8d5a8c1cd3e5a

void test_metis(rtt_dsxx::UnitTest &ut) {
  idx_t nVertices = 10;
  idx_t nWeights = 1;
  idx_t nParts = 2;

  idx_t objval(0);
  std::vector<idx_t> part(nVertices, 0);

  // here's the mesh, there is only one valid cut so the expected result (or a
  // mirror of it) should alays be obtained
  //
  //  0 \       / 6
  //  1 \       / 7
  //  2 - 4 - 5 - 8
  //  3 /       \ 9

  // Indexes of starting points in adjacent array
  idx_t xadj[] = {0, 1, 2, 3, 4, 9, 14, 15, 16, 17, 18};

  // Adjacent vertices in consecutive index order
  // conn. for:     0, 1, 2, 3, 4            , 5,            ,6, 7, 8, 9
  // index:         0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 14,15,16,17
  idx_t adjncy[] = {4, 4, 4, 4, 0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 5, 5, 5, 5};

  // Weights of vertices
  // if all weights are equal then can be set to NULL
  std::vector<idx_t> vwgt(nVertices * nWeights, 0);

  // Partition a graph into k parts using either multilevel recursive
  // bisection or multilevel k-way partitioning.
  int ret =
      METIS_PartGraphKway(&nVertices, &nWeights, xadj, adjncy, NULL, NULL, NULL,
                          &nParts, NULL, NULL, NULL, &objval, &part[0]);

  std::cout << "partition: ";
  for (int32_t i = 0; i < nVertices; ++i) {
    std::cout << part[i] << " ";
  }
  std::cout << std::endl;

  if (ret == METIS_OK)
    PASSMSG("Successfully called METIS_PartGraphKway().");
  else
    FAILMSG("Call to METIS_PartGraphKway() failed.");

  int expectedResult[] = {1, 1, 1, 1, 1, 0, 0, 0, 0, 0};
  int mirrorExpectedResult[] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
  std::vector<idx_t> vExpectedResult(expectedResult,
                                     expectedResult + nVertices);
  std::vector<idx_t> vMirrorExpectedResult(mirrorExpectedResult,
                                           mirrorExpectedResult + nVertices);
  if (part == vExpectedResult || part == vMirrorExpectedResult)
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
