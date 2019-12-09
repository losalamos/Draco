//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/test/tstGatherScatter.cc
 * \author Kent Budge
 * \date   Wed Apr 28 09:31:51 2010
 * \brief  Test c4::gather and c4::scatter functions
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#include "c4/ParallelUnitTest.hh"
#include "c4/gatherv.hh"
#include "c4/scatterv.hh"
#include "ds++/Release.hh"

using namespace std;
using namespace rtt_dsxx;
using namespace rtt_c4;

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstAllgather(UnitTest &ut) {
  int mypid = rtt_c4::node();
  vector<int> allpid(rtt_c4::nodes());
  int status = allgather(&mypid, &allpid[0], 1);
  if (status == 0) {
    PASSMSG("status is okay from allgather");
    status = 0;
    for (int i = 0; i < rtt_c4::nodes(); ++i) {
      if (allpid[i] != i) {
        status = 1;
        break;
      }
    }
    if (status == 0) {
      PASSMSG("gathered values are okay for allgather");
    } else {
      FAILMSG("gathered values are NOT okay for allgather");
    }
  } else
    FAILMSG("status is NOT okay from allgather");
}

//---------------------------------------------------------------------------------------//
void tstDeterminateGatherScatter(UnitTest &ut) {
  unsigned pid = node();
  unsigned const number_of_processors = nodes();
  vector<unsigned> pids(number_of_processors);
  gather(&pid, &pids[0], 1);

  if (pid == 0)
    for (unsigned i = 1; i < number_of_processors; ++i)
      pids[i] += pids[i - 1];

  unsigned base;
  scatter(&pids[0], &base, 1);

  if (base == pid * (pid + 1) / 2)
    PASSMSG("correct base summation in gather/scatter");
  else
    FAILMSG("NOT correct base summation in gather/scatter");

  return;
}

//---------------------------------------------------------------------------//
void tstIndeterminateGatherScatterv(UnitTest &ut) {
  unsigned const pid = node();
  unsigned const number_of_processors = nodes();

  vector<unsigned> send(pid, pid);
  vector<vector<unsigned>> receive;

  indeterminate_gatherv(send, receive);
  PASSMSG("No exception thrown");

  if (pid == 0) {
    if (receive.size() == number_of_processors)
      PASSMSG("correct number of processors in gatherv");
    else
      FAILMSG("NOT correct number of processors in gatherv");
    for (unsigned p = 0; p < number_of_processors; ++p) {
      if (receive[p].size() != p) {
        FAILMSG("NOT correct number of elements in gatherv");
        for (unsigned i = 0; i < p; ++i) {
          if (receive[p][i] != p)
            FAILMSG("NOT correct values in gatherv");
        }
      }
      // Prepare for next test
      receive[p].resize(0);
      receive[p].resize(2 * p, 3 * p);
    }
  }

  indeterminate_scatterv(receive, send);

  if (send.size() == 2 * pid)
    PASSMSG("correct number of processors in scatterv");
  else
    FAILMSG("NOT correct number of processors in scatterv");
  for (unsigned i = 0; i < 2 * pid; ++i) {
    if (send[i] != 3 * pid)
      FAILMSG("NOT correct values in scatterv");
  }

  // Test with empty container
  {
    vector<unsigned> emptysend;
    vector<vector<unsigned>> emptyreceive;

    indeterminate_gatherv(emptysend, emptyreceive);
    PASSMSG(
        "No exception thrown for indeterminate_gatherv with empty containers.");

    if (emptysend.size() != 0)
      ITFAILS;
    if (emptyreceive.size() != number_of_processors)
      ITFAILS;
    if (emptyreceive[pid].size() != 0)
      ITFAILS;

    indeterminate_scatterv(emptyreceive, emptysend);

    if (emptysend.size() != 0)
      ITFAILS;
    if (emptyreceive.size() != number_of_processors)
      ITFAILS;
    if (emptyreceive[pid].size() != 0)
      ITFAILS;
  }

  return;
}

//---------------------------------------------------------------------------//
void tstDeterminateGatherScatterv(UnitTest &ut) {
  unsigned pid = node();
  unsigned const number_of_processors = nodes();
  vector<unsigned> send(pid, pid);
  vector<vector<unsigned>> receive(number_of_processors);
  for (unsigned p = 0; p < number_of_processors; ++p) {
    receive[p].resize(p, p);
  }
  determinate_gatherv(send, receive);

  PASSMSG("No exception thrown");

  if (pid == 0) {
    if (receive.size() == number_of_processors)
      PASSMSG("correct number of processors in gatherv");
    else
      FAILMSG("NOT correct number of processors in gatherv");
    for (unsigned p = 0; p < number_of_processors; ++p) {
      if (receive[p].size() != p) {
        FAILMSG("NOT correct number of elements in gatherv");
        for (unsigned i = 0; i < p; ++i) {
          if (receive[p][i] != p)
            FAILMSG("NOT correct values in gatherv");
        }
      }
      // Prepare for next test
      receive[p].resize(0);
      receive[p].resize(2 * p, 3 * p);
    }
  }

  send.resize(2 * pid);

  determinate_scatterv(receive, send);

  if (send.size() == 2 * pid)
    PASSMSG("correct number of processors in scatterv");
  else
    FAILMSG("NOT correct number of processors in scatterv");
  for (unsigned i = 0; i < 2 * pid; ++i) {
    if (send[i] != 3 * pid)
      FAILMSG("NOT correct values in scatterv");
  }
  return;
}

//---------------------------------------------------------------------------//

void topology_report(UnitTest &ut) {
  size_t const mpi_ranks = rtt_c4::nodes();
  size_t const my_mpi_rank = rtt_c4::node();

  if (my_mpi_rank == 0)
    std::cout << "\nStarting topology_report()..." << std::endl;

  // Store proc name on local proc
  std::string my_pname = rtt_c4::get_processor_name();
  size_t namelen = my_pname.size();

  // Create a container on IO proc to hold names of all nodes.
  vector<std::string> procnames(mpi_ranks);

  // Gather names into pnames on IO proc.
  rtt_c4::indeterminate_gatherv(my_pname, procnames);

  // Look at the data found on the IO proc.
  if (my_mpi_rank == 0) {

    if (procnames[my_mpi_rank].size() != namelen)
      ITFAILS;

    // Count unique processors
    vector<string> unique_processor_names;
    for (size_t i = 0; i < mpi_ranks; ++i) {
      bool found(false);
      for (size_t j = 0; j < unique_processor_names.size(); ++j)
        if (procnames[i] == unique_processor_names[j])
          found = true;
      if (!found)
        unique_processor_names.push_back(procnames[i]);
    }

    // Print a report
    std::cout << "\nWe are using " << mpi_ranks << " mpi rank(s) on "
              << unique_processor_names.size() << " unique nodes.";

    for (size_t i = 0; i < mpi_ranks; ++i) {
      std::cout << "\n  - MPI rank " << i << " is on " << procnames[i];
      if (procnames[i].size() < 1)
        ITFAILS;
    }
    std::cout << std::endl;

    // Generate a map with the node name as the key and a list of MPI
    // ranks as a vector<int> of data.

    vector<vector<size_t>> map_proc_to_ranks(unique_processor_names.size());
    for (size_t i = 0; i < mpi_ranks; ++i) {
      size_t node_number(0);
      for (size_t j = 0; j < unique_processor_names.size(); ++j)
        if (procnames[i] == unique_processor_names[j]) {
          node_number = j;
          break;
        }
      map_proc_to_ranks[node_number].push_back(i);
    }

    std::cout << "\nMPI ranks per node:";
    for (size_t j = 0; j < unique_processor_names.size(); ++j) {
      std::cout << "\n  - Node " << j << " (" << unique_processor_names[j]
                << "): ";
      for (size_t i = 0; i < map_proc_to_ranks[j].size(); ++i)
        std::cout << map_proc_to_ranks[j][i] << ",";
    }
    std::cout << std::endl;
  }
  return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, release);
  try {
    tstAllgather(ut);
    tstDeterminateGatherScatter(ut);
    tstIndeterminateGatherScatterv(ut);
    tstDeterminateGatherScatterv(ut);
    topology_report(ut);
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing tstGatherScatter, " << err.what()
              << endl;
    ut.numFails++;
  } catch (...) {
    std::cout << "ERROR: While testing tstGatherScatter, "
              << "An unknown exception was thrown." << endl;
    ut.numFails++;
  }
  return ut.numFails;
}

//---------------------------------------------------------------------------//
//                        end of tstGatherScatter.cc
//---------------------------------------------------------------------------//
