//----------------------------------*-C++-*----------------------------------//
// config.hh
// Geoffrey Furnish
// Tue Jan 17 10:23:47 1995
//---------------------------------------------------------------------------//
// @> C4 master configuration header file.
//---------------------------------------------------------------------------//

#ifndef __c4_config_hh__
#define __c4_config_hh__

// Start by assuming we are running on a uniprocessor.

#define __C4_SCALAR__

// Now start checking for supported machines.  

// Currently, we assume if it's a Paragon, then they want to use NX.  Running
// MPI on a Paragon would require updating this logic.

#ifdef __PARAGON__
#undef __C4_SCALAR__
#include <nx.h>
#define PGN(a) a
#else
#define PGN(a)
#endif

// Handle MPI ...

#ifdef __LAM__
#define __MPI__
#endif

#ifdef _CRAYMPP
#ifdef C4_SHMEM
#undef __C4_SCALAR__
#include <intrinsics.h>
#include <mpp/shmem.h>
#else
#define __MPI__
#endif
#endif

#if defined(_POWER) && !defined(NOMPI)
#define __MPI__
#endif

#if defined(__sgi) && !defined(NOMPI)
#define __MPI__
#endif

#ifdef __MPI__
#undef __C4_SCALAR__
#ifdef _POWER
extern "C" {
#endif
#include <mpi.h>
#ifdef _POWER
}
#endif
#define MPI(a) a
#else
#define MPI(a)
#endif

// Add one for SCALAR for completeness

#ifdef __C4_SCALAR__
#define SCALAR(a) a
#else
#define SCALAR(a)
#endif

#endif                          // __c4_config_hh__

//---------------------------------------------------------------------------//
//                              end of c4/config.hh
//---------------------------------------------------------------------------//
