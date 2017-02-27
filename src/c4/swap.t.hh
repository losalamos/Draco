//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/swap.t.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:17 2002
 * \brief  C4 MPI template implementation.
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef c4_swap_t_hh
#define c4_swap_t_hh

#include "C4_Functions.hh"
#include "swap.hh"
#include "c4/config.h"
#include "ds++/Assert.hh"

namespace rtt_c4 {
using std::vector;

//---------------------------------------------------------------------------//
// EXCHANGE
//---------------------------------------------------------------------------//

#ifdef C4_MPI

template <class T>
DLL_PUBLIC_c4 void determinate_swap(vector<unsigned> const &outgoing_pid,
                                    vector<vector<T>> const &outgoing_data,
                                    vector<unsigned> const &incoming_pid,
                                    vector<vector<T>> &incoming_data, int tag) {
  Require(outgoing_pid.size() == outgoing_data.size());
  Require(incoming_pid.size() == incoming_data.size());
  Require(&outgoing_data != &incoming_data);

  unsigned incoming_processor_count = incoming_pid.size();
  unsigned outgoing_processor_count = outgoing_pid.size();

  // This block is a no-op for with-c4=scalar.
  // Dito when the vectors are of zero-length.
  {

    // Post the asynchronous sends.
    vector<C4_Req> outgoing_C4_Req(outgoing_processor_count);
    for (unsigned p = 0; p < outgoing_processor_count; ++p) {
      outgoing_C4_Req[p] = rtt_c4::send_async(
          (outgoing_data[p].size() > 0 ? &outgoing_data[p][0] : NULL),
          outgoing_data[p].size(), outgoing_pid[p], tag);
    }

    // Post the asynchronous receives
    vector<C4_Req> incoming_C4_Req(incoming_processor_count);
    for (unsigned p = 0; p < incoming_processor_count; ++p) {
      incoming_C4_Req[p] = receive_async(
          (incoming_data[p].size() > 0 ? &incoming_data[p][0] : NULL),
          incoming_data[p].size(), incoming_pid[p], tag);
    }

    // Wait for all the receives to complete.
    wait_all(incoming_processor_count,
             incoming_processor_count > 0 ? &incoming_C4_Req[0] : NULL);

    // Wait until all the posted sends have completed.
    wait_all(outgoing_processor_count,
             outgoing_processor_count > 0 ? &outgoing_C4_Req[0] : NULL);
  }
  return;
}

//---------------------------------------------------------------------------//
template <class T>
DLL_PUBLIC_c4 void determinate_swap(vector<vector<T>> const &outgoing_data,
                                    vector<vector<T>> &incoming_data, int tag) {
  Require(static_cast<int>(outgoing_data.size()) == rtt_c4::nodes());
  Require(static_cast<int>(incoming_data.size()) == rtt_c4::nodes());
  Require(&outgoing_data != &incoming_data);

  { // This block is a no-op for with-c4=scalar

    unsigned const N = rtt_c4::nodes();

    // Post the asynchronous sends.
    vector<C4_Req> outgoing_C4_Req(N);
    for (unsigned p = 0; p < N; ++p) {
      if (outgoing_data[p].size() > 0) {
        outgoing_C4_Req[p] = rtt_c4::send_async(
            (outgoing_data[p].size() > 0 ? &outgoing_data[p][0] : NULL),
            outgoing_data[p].size(), p, tag);
      }
    }

    // Post the asynchronous receives
    vector<C4_Req> incoming_C4_Req(N);
    for (unsigned p = 0; p < N; ++p) {
      if (incoming_data[p].size() > 0) {
        incoming_C4_Req[p] = receive_async(
            (incoming_data[p].size() > 0 ? &incoming_data[p][0] : NULL),
            incoming_data[p].size(), p, tag);
      }
    }

    // Wait for all the receives to complete.

    wait_all(N, &incoming_C4_Req[0]);

    // Wait until all the posted sends have completed.

    wait_all(N, &outgoing_C4_Req[0]);
  }

  return;
}
//---------------------------------------------------------------------------//
template <class T>
DLL_PUBLIC_c4 void semideterminate_swap(vector<unsigned> const &outgoing_pid,
                                        vector<vector<T>> const &outgoing_data,
                                        vector<unsigned> const &incoming_pid,
                                        vector<vector<T>> &incoming_data,
                                        int tag) {
  Require(outgoing_pid.size() == outgoing_data.size());
  Require(&outgoing_data != &incoming_data);

  unsigned incoming_processor_count = incoming_pid.size();
  unsigned outgoing_processor_count = outgoing_pid.size();

  { // This block is a no-op for with-c4=scalar

    // Send the sizing information using determinate_swap
    vector<vector<unsigned>> outgoing_size(outgoing_pid.size(),
                                           vector<unsigned>(1));
    vector<vector<unsigned>> incoming_size(incoming_pid.size(),
                                           vector<unsigned>(1));

    for (unsigned p = 0; p < outgoing_processor_count; ++p) {
      outgoing_size[p][0] = outgoing_data[p].size();
    }

    determinate_swap(outgoing_pid, outgoing_size, incoming_pid, incoming_size,
                     tag);

    // Post the asynchronous sends.
    vector<C4_Req> outgoing_C4_Req(outgoing_processor_count);
    for (unsigned p = 0; p < outgoing_processor_count; ++p) {
      outgoing_C4_Req[p] = rtt_c4::send_async(
          (outgoing_data[p].size() > 0 ? &outgoing_data[p][0] : NULL),
          outgoing_data[p].size(), outgoing_pid[p], tag);
    }

    // Post the asynchronous receives
    incoming_data.resize(incoming_pid.size());
    vector<C4_Req> incoming_C4_Req(incoming_processor_count);
    for (unsigned p = 0; p < incoming_processor_count; ++p) {
      incoming_data[p].resize(incoming_size[p][0]);
      incoming_C4_Req[p] = receive_async(
          (incoming_data[p].size() > 0 ? &incoming_data[p][0] : NULL),
          incoming_data[p].size(), incoming_pid[p], tag);
    }

    // Wait for all the receives to complete.
    if (incoming_C4_Req.size() > 0)
      wait_all(incoming_processor_count, &incoming_C4_Req[0]);

    // Wait until all the posted sends have completed.
    if (outgoing_C4_Req.size() > 0)
      wait_all(outgoing_processor_count, &outgoing_C4_Req[0]);
  }

  return;
}

//---------------------------------------------------------------------------//
// These functions do nothing if there is no communicator (C4_SCALAR=1)
//---------------------------------------------------------------------------//
#else

template <class T>
DLL_PUBLIC_c4 void determinate_swap(vector<unsigned> const & /*outgoing_pid*/,
                                    vector<vector<T>> const & /*outgoing_data*/,
                                    vector<unsigned> const & /*incoming_pid*/,
                                    vector<vector<T>> & /*incoming_data*/,
                                    int /*tag*/) {
  return;
}
template <class T>
DLL_PUBLIC_c4 void determinate_swap(vector<vector<T>> const & /*outgoing_data*/,
                                    vector<vector<T>> & /*incoming_data*/,
                                    int /*tag*/) {
  return;
}
template <class T>
DLL_PUBLIC_c4 void
semideterminate_swap(vector<unsigned> const & /*outgoing_pid*/,
                     vector<vector<T>> const & /*outgoing_data*/,
                     vector<unsigned> const & /*incoming_pid*/,
                     vector<vector<T>> & /*incoming_data*/, int /*tag*/) {
  return;
}
#endif // C4_MPI

} // end namespace rtt_c4

#endif // c4_swap_t_hh

//---------------------------------------------------------------------------//
// end of c4/swap.t.hh
//---------------------------------------------------------------------------//
