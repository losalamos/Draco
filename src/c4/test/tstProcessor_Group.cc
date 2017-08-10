//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstProcessor_Group.cc
 * \author Kent Budge
 * \date   Tue Sep 21 11:45:44 2004
 * \brief  Unit tests for the Processor_Group class and member functions.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/Processor_Group.hh"
#include "ds++/Release.hh"
#include "ds++/Soft_Equivalence.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

#ifdef C4_MPI
void tstProcessor_Group(rtt_dsxx::UnitTest &ut) {
  unsigned const pid = rtt_c4::node();

  // Test construction
  Processor_Group comm(2);
  PASSMSG("Processor_Group constructed without exception.");

  // Test sum
  unsigned const group_pids = comm.size();
  unsigned const base = pid % 2;
  vector<double> sum(1, base + 1.);
  comm.sum(sum);
  if (rtt_dsxx::soft_equiv(sum[0], group_pids * (base + 1.)))
    PASSMSG("Correct processor group sum");
  else
    FAILMSG("NOT correct processor group sum");

  // Test assemble_vector
  vector<double> myvec;
  size_t const vlen(5);
  for (size_t i = 0; i < vlen; ++i)
    myvec.push_back(pid * 1000 + i);
  vector<double> globalvec;
  comm.assemble_vector(myvec, globalvec);

  if (globalvec.size() != group_pids * vlen)
    ITFAILS;

  // Check the the first 5 elements
  vector<double> goldglobalvec;
  for (size_t j = 0; j < group_pids; ++j)
    for (size_t i = 0; i < vlen; ++i)
      goldglobalvec.push_back((base + 2 * j) * 1000 + i);

  if (!rtt_dsxx::soft_equiv(goldglobalvec.begin(), goldglobalvec.end(),
                            globalvec.begin(), globalvec.end()))
    ITFAILS;

  return;
}
#endif // C4_MPI

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ParallelUnitTest ut(argc, argv, release);
  try {
#ifdef C4_MPI
    tstProcessor_Group(ut);
#else
    ut.passes("Test inactive for scalar");
#endif //C4_MPI
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstProcessor_Group.cc
//---------------------------------------------------------------------------//
