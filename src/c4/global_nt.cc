//----------------------------------*-C++-*----------------------------------//
// global.cc
// Geoffrey Furnish
// Tue Dec 20 1994
//---------------------------------------------------------------------------//
// @> Improves on the NX API for global operations.
//---------------------------------------------------------------------------//

#include "c4/global.hh"

//---------------------------------------------------------------------------//
// This section is just the compatibility wrappers.  The "real" stuff goes
// into the C4 namespace, and is embodied in the implementation files for
// each platform.

// Uhh, maybe these will all end up being inlined...

//---------------------------------------------------------------------------//

#ifdef C4_NX
#include "c4/global_nx.cc"
#endif
#ifdef C4_MPI
#include "c4/global_mpi.cc"
#endif
#ifdef C4_SCALAR
#include "c4/global_scalar.cc"
#endif
#ifdef C4_SHMEM
#include "c4/global_shmem.cc"
#endif

//---------------------------------------------------------------------------//
//                              end of global.cc
//---------------------------------------------------------------------------//
