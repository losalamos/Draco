//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   VendorChecks/test/tstParmetis.cc
 * \date   Monday, May 16, 2016, 16:30 pm
 * \brief  Attempt to link to libparmetis and run a simple problem.
 * \note   Copyright (C) 2016-2019, Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
#include <parmetis.h>
#include <sstream>
#include <vector>

// Example borrowed from the Parmetis manual.

void test_parmetis(rtt_c4::ParallelUnitTest &ut) {
  using std::cout;
  using std::endl;

  Insist(rtt_c4::nodes() == 3,
         "test_parmetis must be called with 3 MPI ranks exactly.");

  // MPI VARIABLES
  size_t const MPI_PROC_ID = rtt_c4::node();
  size_t const MPI_PROC_TOTAL_NUM = rtt_c4::nodes();

  if (MPI_PROC_ID == 0)
    cout << " Parmetis example from LiberLocus." << '\n';
  cout << "I am Proc " << MPI_PROC_ID << '\n';

  // Needed by parmetis

  // These store the (local) adjacency structure of the graph at each
  // processor.
  std::vector<idx_t> xadj(6);
  std::vector<idx_t> adjncy;
  // These store the weights of the vertices and edges.
  idx_t *vwgt = NULL, *adjwgt = NULL;
  // This is used to indicate if the graph is weighted. (0 == no weights)
  idx_t wgtflag = 0;
  // C-style numbering that starts from 0.
  idx_t numflag = 0;
  // This is used to specify the number of weights that each vertex has. It is
  // also the number of balance constraints that must be satisfied.
  idx_t ncon = 1;
  // This is used to specify the number of sub-domains that are desired. Note
  // that the number of subdomains is independent of the number of processors
  // that call this routine.
  idx_t nparts = 3;
  // An array of size ncon × nparts that is used to specify the fraction of
  // vertex weight that should be distributed to each sub-domain for each
  // balance constraint. If all of the sub-domains are to be of the same size
  // for every vertex weight, then each of the ncon × nparts elements should
  // be set to a value of 1/nparts.
  std::vector<real_t> tpwgts(ncon * nparts, static_cast<real_t>(1.0 / nparts));
  // An array of size ncon that is used to specify the imbalance tolerance for
  // each vertex weight, with 1 being perfect balance and nparts being perfect
  // imbalance. A value of 1.05 for each of the ncon weights is recommended.
  real_t ubvec(static_cast<real_t>(1.05));
  // This is an array of integers that is used to pass additional parameters
  // for the routine.
  std::vector<idx_t> options(4, 0);
  // Upon successful completion, the number of edges that are cut by the
  // partitioning is written to this parameter.
  idx_t edgecut(0);

  MPI_Comm_dup(MPI_COMM_WORLD, &rtt_c4::communicator);

  // This is an array of size equal to the number of locally-stored
  // vertices. Upon successful completion the partition vector of the
  // locally-stored vertices is written to this array.
  Check(MPI_PROC_ID < INT_MAX);
  std::vector<idx_t> part(5, static_cast<idx_t>(MPI_PROC_ID));

  // This array describes how the vertices of the graph are distributed among
  // the processors. Its contents are identical for every processor.
  std::vector<idx_t> vtxdist = {0, 5, 10, 15};

  // Dependent on each processor
  if (MPI_PROC_ID == 0) {
    adjncy.resize(13);

    xadj[0] = 0;
    xadj[1] = 2;
    xadj[2] = 5;
    xadj[3] = 8;
    xadj[4] = 11;
    xadj[5] = 13;

    adjncy[0] = 1;
    adjncy[1] = 5;
    adjncy[2] = 0;
    adjncy[3] = 2;
    adjncy[4] = 6;
    adjncy[5] = 1;
    adjncy[6] = 3;
    adjncy[7] = 7;
    adjncy[8] = 2;
    adjncy[9] = 4;
    adjncy[10] = 8;
    adjncy[11] = 3;
    adjncy[12] = 9;
  } else if (MPI_PROC_ID == 1) {
    adjncy.resize(18);

    xadj[0] = 0;
    xadj[1] = 3;
    xadj[2] = 7;
    xadj[3] = 11;
    xadj[4] = 15;
    xadj[5] = 18;

    adjncy[0] = 0;
    adjncy[1] = 6;
    adjncy[2] = 10;
    adjncy[3] = 1;
    adjncy[4] = 5;
    adjncy[5] = 7;
    adjncy[6] = 11;
    adjncy[7] = 2;
    adjncy[8] = 6;
    adjncy[9] = 8;
    adjncy[10] = 12;
    adjncy[11] = 3;
    adjncy[12] = 7;
    adjncy[13] = 9;
    adjncy[14] = 13;
    adjncy[15] = 4;
    adjncy[16] = 8;
    adjncy[17] = 14;
  } else if (MPI_PROC_ID == 2) {
    adjncy.resize(13);

    xadj[0] = 0;
    xadj[1] = 2;
    xadj[2] = 5;
    xadj[3] = 8;
    xadj[4] = 11;
    xadj[5] = 13;

    adjncy[0] = 5;
    adjncy[1] = 11;
    adjncy[2] = 6;
    adjncy[3] = 10;
    adjncy[4] = 12;
    adjncy[5] = 7;
    adjncy[6] = 11;
    adjncy[7] = 13;
    adjncy[8] = 8;
    adjncy[9] = 12;
    adjncy[10] = 14;
    adjncy[11] = 9;
    adjncy[12] = 13;
  }
  if (MPI_PROC_ID == 0)
    cout << "parmetis initialized." << '\n';

  int result = ParMETIS_V3_PartKway(&vtxdist[0], &xadj[0], &adjncy[0], vwgt,
                                    adjwgt, &wgtflag, &numflag, &ncon, &nparts,
                                    &tpwgts[0], &ubvec, &options[0], &edgecut,
                                    &part[0], &rtt_c4::communicator);

  if (result == METIS_OK) {
    std::ostringstream msg;
    msg << "[" << MPI_PROC_ID
        << "] ParMETIS_V3_AdaptiveRepart did not return an error.";
    PASSMSG(msg.str());
  } else {
    std::ostringstream msg;
    msg << "[" << MPI_PROC_ID
        << "] ParMETIS_V3_AdaptiveRepart returned an error code.";
    FAILMSG(msg.str());
  }

  if (MPI_PROC_ID == 0)
    cout << "parmetis finalized." << endl;

  for (size_t pid = 0; pid < MPI_PROC_TOTAL_NUM; ++pid) {
    rtt_c4::global_barrier();
    cout << MPI_PROC_ID << " edgecut " << edgecut << '\n';
    for (int i = 0; i < 5; i++)
      cout << "[" << MPI_PROC_ID << "] " << part[i] << endl;
  }

  return;
}

//---------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_dsxx::release);
  try {
    test_parmetis(ut);
  }
  UT_EPILOG(ut);
}

//---------------------------------------------------------------------------//
// end of tstParmetis.cc
//---------------------------------------------------------------------------//
