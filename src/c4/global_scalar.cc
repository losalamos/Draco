//----------------------------------*-C++-*----------------------------------//
// global_scalar.cc
// Maurice LeBrun
// Wed Feb  1 16:01:58 1995
//---------------------------------------------------------------------------//
// @> Global C4 functions for a scalar architecture
//---------------------------------------------------------------------------//

#include "config.hh"

#include <time.h>

//---------------------------------------------------------------------------//
// Miscellaneous

C4_NAMESPACE_BEG

void Init( int& argc, char **& argv ) {}

void Finalize() {}

struct timespec tsclock;

double Wtime()
{
    int clock =
#if defined(CLOCK_SGI_CYCLE)
        CLOCK_SGI_CYCLE;
#else
    CLOCK_REALTIME;
#endif

    clock_gettime( clock, &tsclock );
    return tsclock.tv_sec + tsclock.tv_nsec*1.0e-9;
}

double Wtick() 
{
    int clock =
#if defined(CLOCK_SGI_CYCLE)
        CLOCK_SGI_CYCLE;
#else
    CLOCK_REALTIME;
#endif

    clock_getres( clock, &tsclock );
    return tsclock.tv_sec + tsclock.tv_nsec*1.0e-9;
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
