//----------------------------------*-C++-*----------------------------------//
// config.hh
// Geoffrey Furnish
// Tue Jan 17 10:23:47 1995
//---------------------------------------------------------------------------//
// @> C4 master configuration header file.
//---------------------------------------------------------------------------//

#ifndef __c4_config_hh__
#define __c4_config_hh__

// Configure some namespace support

#define C4_NAMESPACE_BEG namespace C4 {
#define C4_NAMESPACE_END }

// The basic thing we need to determine here, is what the underlying message
// transport system is going to be.  The following code is going to attempt
// to unravel this.  After we figure out which transport layer is being used, 
// we can go include the right files.

#if !defined(C4_SCALAR) && !defined(C4_MPI ) && \
    !defined(C4_NX) && !defined(C4_SHMEM)

// Okay, user did not say explicitly what to do, try to figure it out by
// infering the properties of the system we're on.

// Start by assuming we are running on a uniprocessor.

#define C4_SCALAR

// Now start checking for supported/auto-recognized machines.  

// Currently, we assume if it's a Paragon, then they want to use NX.  Running
// MPI on a Paragon would require updating this logic.

#ifdef __PARAGON__
#undef C4_SCALAR
#define C4_NX
#endif

// Handle MPI ...

#ifdef __LAM__
#undef C4_SCALAR
#define  C4_MPI
#endif

#endif // !defined(C4_XXX)

//---------------------------------------------------------------------------//
// By this point, the transport layer in use is supposed to be determined.
// Now we define some convenience macros.  Also include any platform
// dependent headers which may be required.
//---------------------------------------------------------------------------//

// Intel's NX

#ifdef C4_NX
#define PGN(a) a
#else
#define PGN(a)
#endif

// MPI

#ifdef C4_MPI
#include <mpi.h>
#define MPI(a) a
#else
#define MPI(a)
#endif

// SHMEM

#ifdef C4_SHMEM
#include "c4/shmem.hh"
#else
#endif

// Add one for SCALAR for completeness

#ifdef C4_SCALAR
#define SCALAR(a) a
#else
#define SCALAR(a)
#endif

// We need support for high resolution clocks.  This seems to be very
// platform dependent, so we use macros to parameterize the options. 

#ifdef __linux
#define C4_USE_GETTIMEOFDAY
#endif

#ifdef __sgi
#define C4_USE_CLOCK_GETTIME
#endif

#endif                          // __c4_config_hh__

//---------------------------------------------------------------------------//
//                              end of c4/config.hh
//---------------------------------------------------------------------------//
