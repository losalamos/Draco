//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/bin/xthi.cc
 * \author Mike Berry <mrberry@lanl.gov>, Kelly Thompson <kgt@lanl.gov>
 * \date   Wednesday, Aug 09, 2017, 11:45 am
 * \brief  Print MPI rank, thread number and core affinity bindings.
 * \note   Copyright (C) 2017 Los Alamos National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

// #define _GNU_SOURCE
// #include <cstdio>
// #include <cstring>
// #include <sched.h>
// #include <unistd.h>

#include "c4_omp.h"
#include "c4/C4_Functions.hh"
#include "ds++/SystemCall.hh"
#include <iomanip>
#include <iostream>
#include <sstream>

namespace rtt_c4 {

//----------------------------------------------------------------------------//
/* Borrowed from util-linux-2.13-pre7/schedutils/taskset.c */
std::string cpuset_to_string(void) {

  // return value;
  std::ostringstream cpuset;

  // local storage
  cpu_set_t coremask;
  (void)sched_getaffinity(0, sizeof(coremask), &coremask);

  size_t entry_made = 0;
  for (int i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, &coremask)) {
      int run = 0;
      entry_made = 1;
      for (int j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, &coremask))
          run++;
        else
          break;
      }
      if (run == 0)
        cpuset << i << ",";
      else
        cpuset << i << "" << i + run << ",";
      i += run;
    }
  }
  return cpuset.str().substr(0, cpuset.str().length() - entry_made);
}

} // end namespace rtt_c4

//----------------------------------------------------------------------------//
int main(int argc, char *argv[]) {

  rtt_c4::initialize(argc, argv);
  int const rank = rtt_c4::node();
  std::string const hostname = rtt_dsxx::draco_gethostname();

#pragma omp parallel
  {
    int thread = omp_get_thread_num();
    std::string cpuset = rtt_c4::cpuset_to_string();

#pragma omp critical
    {
      std::cout << hostname << " :: Rank " << std::setfill('0') << std::setw(5)
                << rank << ", Thread " << std::setfill('0') << std::setw(3)
                << thread << ", core affinity = " << cpuset << std::endl;
    } // end omp critical
  }   // end omp parallel

  rtt_c4::finalize();
  return (0);
}

//----------------------------------------------------------------------------//
// End c4/bin/xthi.cc
//----------------------------------------------------------------------------//
