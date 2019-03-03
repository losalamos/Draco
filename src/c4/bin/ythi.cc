//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/bin/xthi.cc
 * \author Mike Berry <mrberry@lanl.gov>, Kelly Thompson <kgt@lanl.gov>,
 *         Tim Kelley <tkelley@lanl.gov.
 * \date   Tuesday, Jun 05, 2018, 17:12 pm
 * \brief  Print MPI rank, thread number and core affinity bindings.
 * \note   Copyright (C) 2018-2019 Triad National Security, LLC.
 *         All rights reserved.
 *
 * Rewritten by Tim Kelley to run C++11 std::threads You may override
 * \c NUM_WORKERS on the compile command line.  For example to run with 4 worker
 * threads:
 *
 * \code
 * $ ./ythi 4
 * \endcode
 *
 * The default is 1 worker thread (over and above the host thread)
 */
//---------------------------------------------------------------------------//

#include "c4/C4_Functions.hh"
#include "c4/bin/xthi.hh"
#include <atomic>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

//----------------------------------------------------------------------------//
/**\brief After atomic bool changes to true, print out some thread info. */
void run_thread(std::atomic<bool> &signal, std::string const &hostname,
                int const rank, size_t const simple_thread_id) {
  while (!signal) {
  }
  unsigned const num_cpu = std::thread::hardware_concurrency();
  std::string cpuset = rtt_c4::cpuset_to_string(num_cpu);
  std::cout << hostname << " :: Rank " << std::setfill('0') << std::setw(5)
            << rank << ", Thread " << std::setfill('0') << std::setw(3)
            << simple_thread_id << ", core affinity = " << cpuset << std::endl;
  return;
}

//----------------------------------------------------------------------------//
int main(int argc, char **argv) {
  size_t const YTHI_NUM_WORKERS = (argc > 1) ? std::stoi(argv[1]) : 1;
  unsigned const num_cpus = std::thread::hardware_concurrency();

  rtt_c4::initialize(argc, argv);
  int const rank = rtt_c4::node();
  if (rank == 0)
    std::cout << "Found " << num_cpus << " logical CPUs per node." << std::endl;

  std::string const hostname = rtt_dsxx::draco_gethostname();

  std::vector<std::atomic<bool>> signals(YTHI_NUM_WORKERS);
  std::vector<std::thread> threads(YTHI_NUM_WORKERS);

  for (size_t i = 0; i < YTHI_NUM_WORKERS; ++i) {
    signals[i].store(false);
    threads[i] = std::thread(run_thread, std::ref(signals[i]),
                             std::ref(hostname), rank, i + 1);
  }
  std::string cpuset = rtt_c4::cpuset_to_string(num_cpus);
  int const host_thread(0);
  std::cout << hostname << " :: Rank " << std::setfill('0') << std::setw(5)
            << rank << ", Thread " << std::setfill('0') << std::setw(3)
            << host_thread << ", core affinity = " << cpuset << std::endl;
  for (size_t i = 0; i < YTHI_NUM_WORKERS; ++i) {
    signals[i].store(true);
    threads[i].join();
  }

  rtt_c4::finalize();
  return (0);
}

//----------------------------------------------------------------------------//
// End c4/bin/ythi.cc
//----------------------------------------------------------------------------//
