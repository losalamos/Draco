//----------------------------------*-C++-*----------------------------------//
// global_scalar.cc
// Maurice LeBrun
// Wed Feb  1 16:01:58 1995
//---------------------------------------------------------------------------//
// @> Global C4 functions for a scalar architecture
//---------------------------------------------------------------------------//

#include "config.hh"

#include <time.h>

#ifdef C4_USE_GETTIMEOFDAY
#include <sys/time.h>
#include <unistd.h>
#endif

//---------------------------------------------------------------------------//
// Miscellaneous

C4_NAMESPACE_BEG

void Init( int& argc, char **& argv ) {}

void Finalize() {}

#ifdef C4_USE_CLOCK_GETTIME
struct timespec tsclock;
#endif
#ifdef C4_USE_GETTIMEOFDAY
struct timeval tvclock;
#endif


double Wtime()
{
#ifdef C4_USE_CLOCK_GETTIME
    int clock =
#if defined(CLOCK_SGI_CYCLE)
        CLOCK_SGI_CYCLE;
#else
    CLOCK_REALTIME;
#endif

    clock_gettime( clock, &tsclock );
    return tsclock.tv_sec + tsclock.tv_nsec*1.0e-9;
#endif
#ifdef C4_USE_GETTIMEOFDAY
    gettimeofday( &tvclock, NULL );
    return tvclock.tv_sec + tvclock.tv_usec*1.0e-6;
#endif
}

double Wtick() 
{
#ifdef C4_USE_CLOCK_GETTIME
    int clock =
#if defined(CLOCK_SGI_CYCLE)
        CLOCK_SGI_CYCLE;
#else
    CLOCK_REALTIME;
#endif

    clock_getres( clock, &tsclock );
    return tsclock.tv_sec + tsclock.tv_nsec*1.0e-9;
#endif
#ifdef C4_USE_GETTIMEOFDAY
    return 1.0e-3;              // Oh, this is so bogus!
#endif
}

int  node()
{
    return 0;
}

int  nodes()
{
    return 1;
}

int  group()
{
    return 0;
}

void gsync() {}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of global_scalar.cc
//---------------------------------------------------------------------------//
