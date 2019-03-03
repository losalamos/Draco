//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_MPI.hh
 * \author Thomas M. Evans
 * \date   Thu Mar 21 16:56:16 2002
 * \brief  C4 MPI function declarations.
 * \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#ifndef c4_C4_MPI_hh
#define c4_C4_MPI_hh

#include "c4/config.h"
#include <algorithm>

#ifdef C4_MPI

#include "MPI_Traits.hh"
#include "c4_mpi.h"
#include "ds++/Assert.hh"

namespace rtt_c4 {

//---------------------------------------------------------------------------//
// MPI Communicator
//---------------------------------------------------------------------------//

DLL_PUBLIC_c4 extern MPI_Comm communicator;
extern bool initialized;

//---------------------------------------------------------------------------//
// SETUP FUNCTIONS
//---------------------------------------------------------------------------//

template <typename Comm> void inherit(const Comm &comm) {
  Remember(int result =) MPI_Comm_dup(comm, &communicator);
  Check(result == MPI_SUCCESS);
}

template <class T>
int create_vector_type(unsigned count, unsigned blocklength, unsigned stride,
                       C4_Datatype &new_type) {
  int info = MPI_Type_vector(count, blocklength, stride,
                             MPI_Traits<T>::element_type(), &new_type);

  if (info != C4_SUCCESS)
    return info;

  info = MPI_Type_commit(&new_type);

  return info;
}

//----------------------------------------------------------------------------//
/*!
 * Broadcast the range [first, last) from proc 0 into [result, ...) on all other
 * processors.
 */

// This signature must be exported since it is explicitly instantiated.
template <typename T>
DLL_PUBLIC_c4 int broadcast(T * /*buffer*/, int /*size*/, int /*root*/);

// This signature is defined in the header so no export is required.
template <typename ForwardIterator, typename OutputIterator>
void broadcast(ForwardIterator first, ForwardIterator last,
               OutputIterator result) {
  typedef typename std::iterator_traits<ForwardIterator>::value_type value_type;
  typedef
      typename std::iterator_traits<ForwardIterator>::difference_type diff_type;

  // Proc 0 does not copy any data into the result iterator.

  diff_type size;
  if (node() == 0)
    size = std::distance(first, last);

  Remember(int check =) broadcast(&size, 1, 0);
  Check(check == MPI_SUCCESS);

  value_type *buf = new value_type[size];
  if (node() == 0)
    std::copy(first, last, buf);

  Check(size < INT_MAX);
  Remember(check =) broadcast(buf, static_cast<int>(size), 0);
  Check(check == MPI_SUCCESS);

  if (node() != 0)
    std::copy(buf, buf + size, result);

  delete[] buf;
  return;
}

// safer version of broadcast using stl ranges
template <typename ForwardIterator, typename OutputIterator>
void broadcast(ForwardIterator first, ForwardIterator last,
               OutputIterator result, OutputIterator result_end) {
  // Check that the result is large enough to hold the data that is currently in
  // buf.
  Insist(std::distance(first, last) == std::distance(result, result_end),
         "Destination must be same size as source data.");
  broadcast(first, last, result);
  return;
}

} // end namespace rtt_c4

#include "C4_MPI.i.hh"

#endif // C4_MPI

#endif // c4_C4_MPI_hh

//---------------------------------------------------------------------------//
// end of c4/C4_MPI.hh
//---------------------------------------------------------------------------//
