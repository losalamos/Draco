//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Invert_Comm_Map.cc
 * \author Mike Buksas, Rob Lowrie
 * \date   Mon Nov 19 10:09:11 2007
 * \brief  Implementation of Invert_Comm_Map
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//

#include "Invert_Comm_Map.hh"
#include "MPI_Traits.hh"
#include "ds++/Assert.hh"
#include <vector>

namespace rtt_c4 {

//---------------------------------------------------------------------------//
// MPI version of invert_comm_map
#ifdef C4_MPI
void invert_comm_map(Invert_Comm_Map_t const &to_map,
                     Invert_Comm_Map_t &from_map) {
  const int myproc = rtt_c4::node();
  const int numprocs = rtt_c4::nodes();

  // The local vector that the other procs will set the size they are sending,
  // at the index of their processor number.  Zero indicates no comm from that
  // proc.
  std::vector<size_t> proc_flag(numprocs, 0); // initially, all zero

  // Create the RMA memory window of the vector.
  MPI_Win win;
  MPI_Win_create(&proc_flag[0], numprocs * sizeof(size_t), sizeof(size_t),
                 MPI_INFO_NULL, MPI_COMM_WORLD, &win);

  // Assertion value for fences.  Currently, we effectively don't set
  // anything (zero).
  const int fence_assert = 0;

  // Set the local and remote vector values
  MPI_Win_fence(fence_assert, win);
  for (auto it = to_map.begin(); it != to_map.end(); ++it) {
    Require(it->first >= 0);
    Require(it->first < numprocs);
    Require(it->second > 0);
    if (it->first == myproc) {
      // ... set our local value
      proc_flag[myproc] = it->second;
    } else {
      // ... set the value on the remote proc
      MPI_Put(&(it->second), 1, MPI_Traits<size_t>::element_type(), it->first,
              myproc, 1, MPI_Traits<size_t>::element_type(), win);
    }
  }
  MPI_Win_fence(fence_assert, win);

  // Back out the map from the vector
  from_map.clear();
  for (int i = 0; i < numprocs; ++i) {
    if (proc_flag[i] > 0)
      from_map[i] = proc_flag[i];
  }

  MPI_Win_free(&win);
  return;
}
//---------------------------------------------------------------------------//
// SCALAR version of invert_comm_map
#elif defined(C4_SCALAR)
void invert_comm_map(Invert_Comm_Map_t const &to_map,
                     Invert_Comm_Map_t &from_map) {
  Require(to_map.size() == 0u || (to_map.size() == 1u && to_map.at(0) > 0));
  from_map.clear();
  auto it = to_map.find(0);
  if (it != to_map.end()) {
    from_map[0] = it->second;
  }
}
#else
//---------------------------------------------------------------------------//
// Default version of invert_comm_map, which throws an error.
void invert_comm_map(Invert_Comm_Map_t const &, Invert_Comm_Map_t &) {
  Insist(0, "invert_comm_map not implemented for this communication type!");
}
#endif // ifdef C4_MPI

} // end namespace rtt_c4

//---------------------------------------------------------------------------//
// end of Invert_Comm_Map.cc
//---------------------------------------------------------------------------//
