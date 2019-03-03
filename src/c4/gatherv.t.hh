//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/gatherv.t.hh
 * \author Thomas M. Evans, Kent Budge
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI template implementation.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef c4_gatherv_t_hh
#define c4_gatherv_t_hh

#include "C4_Functions.hh"
#include "gatherv.hh"
#include "c4/config.h"
#include "ds++/Assert.hh"
#include <algorithm>
#include <limits>

namespace rtt_c4 {

//---------------------------------------------------------------------------//
// GATHER
//---------------------------------------------------------------------------//

void indeterminate_gatherv(std::string &outgoing_data,
                           std::vector<std::string> &incoming_data) {
  // convert from string to vector<char>
  std::vector<char> outgoing_data_vc(outgoing_data.size());
  std::copy(outgoing_data.begin(), outgoing_data.end(),
            outgoing_data_vc.begin());

  // create container for return value
  std::vector<std::vector<char>> incoming_data_vvc(nodes());

  // Call intederminate_gatherv(...) for vector<char>.
  indeterminate_gatherv(outgoing_data_vc, incoming_data_vvc);

  if (node() == 0) {
    // Move data from vector<vector<char> > back to vector<string> format
    incoming_data.resize(nodes());
    for (size_t i = 0; i < static_cast<size_t>(nodes()); ++i) {
      size_t const len = incoming_data_vvc[i].size();
      incoming_data[i].resize(len);
      std::copy(incoming_data_vvc[i].begin(), incoming_data_vvc[i].end(),
                incoming_data[i].begin());
      Check(incoming_data[i].size() > 0);
    }
  }

  return;
}

//---------------------------------------------------------------------------//
template <typename T>
void indeterminate_gatherv(std::vector<T> &outgoing_data,
                           std::vector<std::vector<T>> &incoming_data) {
#ifdef C4_MPI
  { // This block is a no-op for with-c4=scalar

    unsigned const N(rtt_c4::nodes());
    incoming_data.resize(N);

    Check(outgoing_data.size() < INT_MAX);
    int count(static_cast<int>(outgoing_data.size()));
    if (rtt_c4::node() == 0) {
      std::vector<int> counts(N, -1);
      std::vector<int> displs(N, -1);

      // for valid comm world, this should always be true.
      Check(counts.size() > 0);
      Remember(int check =) gather(&count, &counts[0], 1);
      Check(check == MPI_SUCCESS);
      uint64_t total_count_64(0);
      for (unsigned p = 0; p < N; ++p) {
        Check(total_count_64 < INT_MAX);
        displs[p] = static_cast<int>(total_count_64);
        total_count_64 += counts[p];
      }
      // Require that total_count_64 can be expressed as a 32-bit integer.
      Insist(total_count_64 < std::numeric_limits<int>::max(),
             "The size of the array (nranks*vector.size()) exceeds "
             "32-bit unsigned integer limit.");
      unsigned total_count = static_cast<unsigned>(total_count_64);

      // We can only use the &vec[0] notation if the vector is non-zero
      // in length.  An shorthand-if is used to pass 'NULL' to mpi if
      // there is no data to gather.
      std::vector<T> recbuf(total_count, 42);
      Check(outgoing_data.size() < INT_MAX);
      Remember(check =) rtt_c4::gatherv((count > 0 ? &outgoing_data[0] : NULL),
                                        static_cast<int>(outgoing_data.size()),
                                        (total_count > 0 ? &recbuf[0] : NULL),
                                        &counts[0], &displs[0]);
      Check(check == MPI_SUCCESS);

      for (unsigned p = 0; p < N; ++p) {
        incoming_data[p].assign(recbuf.begin() + displs[p],
                                recbuf.begin() + displs[p] + counts[p]);
      }

    } else {
      Remember(int check =) gather(&count, static_cast<int *>(NULL), 1);
      Check(check == MPI_SUCCESS);
      Remember(check =)
          gatherv((count > 0 ? &outgoing_data[0] : NULL), count,
                  static_cast<T *>(NULL), 0, static_cast<int *>(NULL));
      Check(check == MPI_SUCCESS);
    }
  }
#else
  {
    // Only need to copy outgoing to incoming
    incoming_data.resize(0);
    incoming_data.resize(1, outgoing_data);
  }
#endif // C4_MPI

  return;
}

//---------------------------------------------------------------------------//
template <typename T>
void determinate_gatherv(std::vector<T> &outgoing_data,
                         std::vector<std::vector<T>> &incoming_data) {
  Require(static_cast<int>(incoming_data.size()) == rtt_c4::nodes());

#ifdef C4_MPI
  { // This block is a no-op for with-c4=scalar
    unsigned const N(rtt_c4::nodes());

    Check(outgoing_data.size() < INT_MAX);
    int count(static_cast<int>(outgoing_data.size()));
    if (rtt_c4::node() == 0) {
      std::vector<int> counts(N, -1);
      std::vector<int> displs(N, -1);
      uint64_t total_count_64(0);
      for (unsigned p = 0; p < N; ++p) {
        Check(incoming_data[p].size() < INT_MAX);
        counts[p] = static_cast<int>(incoming_data[p].size());
        displs[p] = static_cast<int>(total_count_64);
        total_count_64 += counts[p];
      }
      // Require that total_count_64 can be expressed as a 32-bit integer.
      Insist(total_count_64 < std::numeric_limits<int>::max(),
             "The size of the array (nranks*vector.size()) exceeds "
             "32-bit unsigned integer limit.");
      unsigned total_count = static_cast<unsigned>(total_count_64);

      std::vector<T> recbuf(total_count, 42);
      // &vec[0] is only valid if vector has non-zero length
      rtt_c4::gatherv((count > 0 ? &outgoing_data[0] : NULL), count,
                      (total_count > 0 ? &recbuf[0] : NULL), &counts[0],
                      &displs[0]);

      for (unsigned p = 0; p < N; ++p) {
        incoming_data[p].resize(counts[p]);
        std::copy(recbuf.begin() + displs[p],
                  recbuf.begin() + displs[p] + counts[p],
                  incoming_data[p].begin());
      }

    } else {
      gatherv((count > 0 ? &outgoing_data[0] : NULL), count,
              static_cast<T *>(NULL), 0, static_cast<int *>(NULL));
    }
  }
#else
  {
    // Only need to copy outgoing to incoming
    incoming_data.resize(0);
    incoming_data.resize(1, outgoing_data);
  }
#endif // C4_MPI

  return;
}

} // end namespace rtt_c4

#endif // c4_gatherv_t_hh

//---------------------------------------------------------------------------//
// end of c4/gatherv.t.hh
//---------------------------------------------------------------------------//
