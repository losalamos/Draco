//----------------------------------*-C++-*----------------------------------//
// timer.hh
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Routine for CPU and Wallclock timing.
//---------------------------------------------------------------------------//

#ifndef __sn_timer_hh__
#define __sn_timer_hh__

#include "sn/precision.hh"

#include <sys/types.h>
#include <time.h>

void timer( REAL &cpu, REAL &wall )
{
    clock_t cpu_temp;
    time_t  wall_temp;
    time_t  *tloc;

    tloc = 0;

    cpu_temp  = clock();
    wall_temp = time( tloc );

    cpu  = cpu_temp / 1000000.0;
    wall = wall_temp;
}

#endif                          // __sn_timer_hh__

//---------------------------------------------------------------------------//
//                              end of timer.hh
//---------------------------------------------------------------------------//

